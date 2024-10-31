
#' Feature Selection and Model Training with Priority Methylation Method
#'
#' This function performs feature selection through iterative Boruta phases, followed by UMAP dimensionality reduction,
#' machine learning model training, and heatmap visualization of methylation beta values. The data should have features in columns and samples in rows, with the label as the dependent variable.
#'
#' @param data A data frame where features are in columns and samples are in rows. The dependent variable (label) should be included in the data.
#' @param label_col A character string specifying the name of the column containing the dependent variable (default is "label").
#' @param split_ratio A numeric value for the ratio of data to be used for training (default is 0.7).
#' @param models A character vector specifying machine learning models to be used (default is "rf" for random forest).
#' @param max_runs_phase1 An integer specifying the maximum number of runs for the first phase of Boruta (default is 100).
#' @param max_runs_phase2 An integer specifying the maximum number of runs for the second phase of Boruta (default is 100).
#' @param max_runs_phase3 An integer specifying the maximum number of runs for the third phase of Boruta (default is 100).
#'
#' @return A list containing:
#' \describe{
#'   \item{models}{A list with trained models and their corresponding confusion matrices.}
#'   \item{umap_df}{A data frame with UMAP results for visualization.}
#'   \item{final_features}{A data frame with the final selected features.}
#'   \item{heatmap}{A heatmap object showing methylation beta values.}
#'   \item{umap_plot}{A ggplot object of the UMAP plot for unsupervised clustering.}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Run priority_meth with default parameters
#' result <- priority_meth(data = your_data, label_col = "your_label_column")
#'
#' # Access UMAP plot
#' print(result$umap_plot)
#'
#' # Access model accuracy
#' print(result$models$rf$confusion_matrix)
#' }

priority_meth <- function(data = data, label_col = "label", split_ratio = 0.7, models = c("rf"),
                          max_runs_phase1 = 100, max_runs_phase2 = 100, max_runs_phase3 = 100) {

    library(Boruta)
    library(dplyr)
    library(umap)
    library(ggplot2)
    library(caret)
    library(pheatmap)
    library(tibble)
    library(plyr)

    # Convert label to factor if not already
    data[[label_col]] <- as.factor(data[[label_col]])

    seeds = c(738, 521, 822)

    # Step 1: Data partitioning (Train-test split)
    set.seed(seeds[1])  # Use first seed for data splitting
    partitionrule <- createDataPartition(data[[label_col]], p = split_ratio, list = FALSE)
    training_set <- data[partitionrule, ]
    testing_set <- data[-partitionrule, ]

    # Function to run Boruta with different seeds and combine results
    run_boruta_repeated <- function(data, label_col, seeds, max_runs_phase1, max_runs_phase2, max_runs_phase3) {

        feature_importance_list <- list()

        for (seed in seeds) {
            set.seed(seed)

            #Phase 1
            boruta_run <- Boruta(as.formula(paste(label_col, "~ .")), data = data, doTrace = 3, maxRuns = max_runs_phase1)
            attstats <- attStats(boruta_run)
            attstats <- filter(attstats, decision != "Rejected")
            attstats <- rownames_to_column(attstats, var = "Sites")

            #Phase 2 Prep
            data_phase_1 <- data %>% select(all_of(c(label_col, attstats$Sites)))

            #Phase 2
            boruta_run_2 <- Boruta(as.formula(paste(label_col, "~ .")), data = data_phase_1, doTrace = 3, maxRuns = max_runs_phase2)
            attstats_2 <- attStats(boruta_run_2)
            attstats_2 <- filter(attstats_2, decision != "Rejected")
            attstats_2 <- rownames_to_column(attstats_2, var = "Sites")

            #Phase 3 Prep
            data_phase_2 <- data_phase_1 %>% select(all_of(c(label_col, attstats_2$Sites)))

            #Phase 3
            boruta_run_3 <- Boruta(as.formula(paste(label_col, "~ .")), data = data_phase_2, doTrace = 3, maxRuns = max_runs_phase3)
            attstats_3 <- attStats(boruta_run_3)
            attstats_3 <- filter(attstats_3, decision == "Confirmed")
            attstats_3 <- rownames_to_column(attstats_3, var = "Sites")

            feature_importance_list <- append(feature_importance_list, list(attstats_3))
        }

        # Combine results by selecting features that are confirmed in all runs
        common_features <- join(feature_importance_list[[1]], feature_importance_list[[2]], by = "Sites")
        common_features <- join(common_features, feature_importance_list[[3]], by = "Sites")
        common_features <- na.omit(common_features)

        return(common_features)
    }

    # Phase 3: Repeated Boruta with different seeds
    final_features <- run_boruta_repeated(training_set, label_col, seeds, max_runs_phase1, max_runs_phase2, max_runs_phase3)
    load("data/cpg_annotations.rda")
    final_features <- join(final_features, cpg_annotations, by = "Sites")

    # Apply selected features to both training and testing sets
    training_set_final <- training_set %>% select(all_of(c(label_col, final_features$Sites)))
    testing_set_final <- testing_set %>% select(all_of(c(label_col, final_features$Sites)))

    # Step 3: UMAP Dimensionality Reduction
    umap_fit <- training_set_final %>%
        select(where(is.numeric)) %>%
        scale() %>%
        umap()

    umap_df <- data.frame(UMAP1 = umap_fit$layout[, 1], UMAP2 = umap_fit$layout[, 2], Label = training_set_final[[label_col]])

    umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Label)) +
        geom_point(size = 2) +
        labs(x = "UMAP1", y = "UMAP2", subtitle = "Unsupervised clustering") +
        theme_minimal()

    # Print UMAP plot
    print(umap_plot)

    # Step 4: Machine learning model training and evaluation
    model_results <- list()

    for (model in models) {
        set.seed(seeds[1])  # Use the first seed for training models
        validation_rule <- trainControl(method = "repeatedcv", number = 10, repeats = 30, classProbs = TRUE)

        trained_model <- train(as.formula(paste(label_col, "~ .")), data = training_set_final, method = model, trControl = validation_rule, metric = "Accuracy")

        # Predict on test set
        test_pred <- predict(trained_model, newdata = testing_set_final)

        # Confusion Matrix
        cm <- confusionMatrix(test_pred, testing_set_final[[label_col]])

        # Store results
        model_results[[model]] <- list(model = trained_model, confusion_matrix = cm)

        # Print results
        print(trained_model)
        print(cm)
    }

    # Step 5: Heatmap
    data_ordered <- training_set_final %>% arrange(training_set_final[[label_col]])

    # Transpose the data for the heatmap (excluding the class column)
    data_ordered_transposed <- as.data.frame(t(data_ordered[,-1]))

    # Create annotation based on the class (dependent variable)
    annotation_col <- data.frame(Class = data_ordered[[label_col]])

    # Set row names of the annotation to match the column names of transposed data
    rownames(annotation_col) <- colnames(data_ordered_transposed)

    # Generate the heatmap with class annotations
    heatmap <- pheatmap(as.matrix(data_ordered_transposed),
                        main = "Heatmap of Methylation Beta Values",
                        cluster_rows = FALSE,
                        cluster_cols = FALSE,
                        show_rownames = TRUE,
                        annotation_col = annotation_col)

    # Output results
    return(list(models = model_results, umap_df = umap_df, final_features = final_features, heatmap = heatmap, umap_plot = umap_plot))
}

