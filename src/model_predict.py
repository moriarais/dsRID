import sys
# Add custom utility scripts directory to the system path for module imports
sys.path.insert(1, '/u/home/r/ryo10244/Xiao_lab/dsRNA_pred/script/utils')

# Import necessary libraries
import re
import pysam
import argparse as ap
import pandas as pd
import numpy as np
from utils import *
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, classification_report, roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC

# Define the main function that handles the primary operations
def main(args):
    print("model_predict.py init")
    
    # Set up argument parser for handling command-line arguments
    argp = ap.ArgumentParser(description="extract features for dsrna prediction from randomly sampled region",
                             formatter_class=ap.ArgumentDefaultsHelpFormatter)
    
    # Add arguments for input files, output directories, and model type
    argp.add_argument(
        "-i", "--input_file",
        help="input file of training and validation data",
        type=str,
        default="/data/dsRID_data.tsv"
    )
    argp.add_argument(
        "-p", "--pred_file",
        help="input file of prediction data",
        type=str,
        default="/data/dsRID_whole.tsv"
    )
    argp.add_argument(
        "-o", "--out_dir",
        help="output directory location",
        type=str,
        default="/data/"
    )
    argp.add_argument(
        "-m", "--model",
        help="model type",
        type=str,
        default="randomf"
    )

    # Parse command-line arguments
    args = argp.parse_args(args)

    # Select model type based on user input
    modeltype = args.model
    model_dic = {
        "logireg": LogisticRegression(random_state=0),
        "randomf": RandomForestClassifier(random_state=0, bootstrap=True, 
                                          criterion="gini", max_features=0.2, min_samples_leaf=4,
                                          min_samples_split=2, n_estimators=100),
        "svm": LinearSVC(random_state=0, max_iter=2000)
    }

    # Initialize the model based on user choice
    model = model_dic[modeltype]
    print("modeltype")
    print(modeltype)

    # Load training and validation data, fill missing values with zero
    all_frame = pd.read_csv(args.input_file, sep='\t')
    all_frame = all_frame.fillna(0)

    # Optional code for feature selection based on distance
    #all_frame['Alu'] = all_frame['dist'].apply(lambda x : 1 if x==0 else 0)
    # Sample positive and negative examples for balanced training
    # all_alt_frame = all_frame.loc[all_frame['label'] == 1].sample(35000)
    # all_null_frame = all_frame.loc[all_frame['label'] == 0].sample(35000)

    # Split positive and negative examples into training sets
    # train_alt_frame = all_alt_frame.iloc[:30000]
    # train_null_frame = all_null_frame.iloc[:30000]

    # Split data into training and validation sets
    # X_train, X_test, y_train, y_test
    train_X, val_X, train_y, val_y = train_test_split(all_frame, all_frame['label'], test_size=0.33, random_state=42)

    # Define columns to be used in training, excluding certain metadata columns
    cols = ~train_X.columns.isin(["mean_start", "mean_end", "coverage", "num_skip", "name", "chr", "start", "end", "label"])
    print('train_X.loc[:, cols]')
    print(train_X.loc[:, cols])

    # Train the model with training data
    model = model.fit(train_X.loc[:, cols], train_y)

    # Optional code for training score and predictions
    #train_sc = model.score(train_X, train_y)
    #train_pred = model.predict_proba(train_X.loc[:, cols])
    #val_pred = model.predict_proba(val_X)

    # Print column names to check selected features
    print(train_X.columns[cols])

    # Calculate feature importance using permutation importance and save to DataFrame (2.6)
    feat_imp = permutation_importance(model, val_X.loc[:, cols], val_y, scoring='r2')
    feat_imp_df = pd.DataFrame(data={'feat': train_X.columns[cols], "mean_R2": feat_imp.importances_mean,
                                     'std_R2': feat_imp.importances_std})
    print(feat_imp_df)

    # Save feature importance to an output file
    feat_imp_df.to_csv(args.out_dir + "/feat_importance_{}.tsv".format(modeltype), sep='\t', index=False)

    # Perform cross-validation and save cross-validation scores
    # evaluate the performance of a model by splitting the data into multiple subsets, or "folds," to ensure that the model is tested on different parts of the data. It helps assess how well the model generalizes to unseen data and reduces the risk of overfitting or underfitting.
    scores = cross_val_score(model, train_X.loc[:, cols], train_y, cv=5)
    sc_frame = pd.DataFrame(data={"scores": scores})
    sc_frame.to_csv(args.out_dir + "/cv_scores_{}.tsv".format(modeltype), sep='\t', index=False)
    print("cross-validation scores")
    print(scores)
    mean_score = np.mean(scores)
    print("Average Cross-Validation Score:", mean_score)

    # Generate predictions for the training set and add them as new columns
    train_pred = model.predict_proba(train_X.loc[:, cols])
    train_X['pred_0'] = train_pred[:, 0]
    train_X['pred_1'] = train_pred[:, 1]

    # Generate probability predictions for the validation set
    val_pred_prob = model.predict_proba(val_X.loc[:, cols])[:, 1]  # Probabilities for Class 1

    # Compute the ROC-AUC score
    roc_auc = roc_auc_score(val_y, val_pred_prob)
    print("ROC-AUC Score for the validation set: {:.2f}".format(roc_auc))

    # Generate predictions for the validation set
    val_pred = model.predict(val_X.loc[:, cols])

    # Evaluate and print the classification report
    class_report = classification_report(val_y, val_pred)
    print("Classification Report for the validation set:")
    print(class_report)

    # Load prediction data, fill missing values, and convert column names to lowercase
    pred_frame = pd.read_csv(args.pred_file, sep='\t')
    print("args.pred_file: " + args.pred_file)
    pred_frame = pred_frame.fillna(0)
    pred_frame.columns = map(str.lower, pred_frame.columns)

    # Define columns to use for prediction
    cols = pred_frame.columns.isin(['std_start', 'std_end', 'len_skip', 'skip_ratio', 'group_num',
                                     'group_std', 'gc_skip', 'bp_start_ct', 'bp_end_ac', 'bp_start_tc',
                                     'bp_start_at', 'bp_end_ag', 'bp_end_tt', 'bp_start_gg', 'bp_start_aa',
                                     'bp_start_tt', 'bp_start_ac', 'bp_start_ag', 'bp_end_tg', 'bp_end_gc',
                                     'bp_end_aa', 'bp_end_gg', 'bp_start_cc', 'bp_start_ca', 'bp_start_gt',
                                     'bp_end_gt', 'bp_end_tc', 'bp_end_cc', 'bp_end_cg', 'bp_end_ct',
                                     'bp_start_gc', 'bp_start_cg', 'bp_end_at', 'bp_end_ta', 'bp_start_tg',
                                     'bp_start_ga', 'bp_start_ta', 'bp_end_ca', 'bp_end_ga'])
                                     # 'bp_end_a', 'bp_end_g', 'bp_end_t', 'bp_start_a', 'bp_start_g', 'bp_start_c', 'bp_end_c', 'bp_end_', 'bp_start_', 'bp_start_t'])

    # Generate predictions for the entire dataset and add them as new columns
    whole_pred = model.predict_proba(pred_frame.loc[:, cols])
    pred_frame['pred_0'] = whole_pred[:, 0]
    pred_frame['pred_1'] = whole_pred[:, 1]
    print("model_predict.py complete")
    print("whole_pred")
    print(whole_pred)
    print("pred_frame")
    print(pred_frame)
    
    # Save prediction results to an output file
    pred_frame.to_csv(args.out_dir + "dsRID_whole.tsv", sep='\t', index=False)

# Run the main function if the script is executed directly
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
