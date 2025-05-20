import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split, cross_val_predict, StratifiedKFold, GroupKFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from scipy import stats

# Load or generate some example data
# Replace this with your actual data loading code

lasso_df = pd.read_csv('/Users/katarinawilloch/Desktop/UiO/Project 1/Data/Initial cleansing/all_metrices_all_labs.csv')

for drug in lasso_df['drug'].unique():
    print(drug)
    lasso_df1 = lasso_df[lasso_df['drug'] == drug]
    print(lasso_df)
    X = lasso_df1[['Tissue', 'Disease.status', 'medium',
        'positive_control', 'transport_time', 'time_until_sample_usage',
        'temperature_frozen ', 'centrifugation_procedure', 'buffer',
        'cell_counting_method', 'microenvironmental_stimuli',
        'cell_culturing_conditions', 'cells', 'plate_reader']]
    y = lasso_df1['DSS2']
    groups = lasso_df1['lab']  # Grouping variable

    # Check if all y values are positive for Box-Cox (add a constant if necessary)
    if (y <= 0).any():
        y += abs(min(y)) + 1  # Shift y to make all values positive

    # Apply Box-Cox transformation
    y_boxcox, lam = stats.boxcox(y)

    scaler_y = StandardScaler()
    y_scaled = scaler_y.fit_transform(y_boxcox.reshape(-1, 1)).flatten() 


    # Convert categorical variables into dummy/one-hot encoded variables
    X_dummies = pd.get_dummies(X, drop_first=True)

    # Scale the features for Lasso regression
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_dummies)

    # Create Lasso model
    lasso_model = Lasso(alpha=0.1)  # Adjust alpha as needed for regularization

    # Cross-validation using GroupKFold
    group_kfold = GroupKFold(n_splits=4)

    # Predict with cross-validation using scaled y
    y_pred_cv_scaled = cross_val_predict(lasso_model, X_scaled, y_scaled, cv=group_kfold, groups=groups)

    # Inverse scaling for y predictions: rescale to the original Box-Cox scale per group
    y_pred_boxcox = scaler_y.inverse_transform(y_pred_cv_scaled.reshape(-1, 1)).flatten()

    # Manually apply the inverse Box-Cox transformation to predictions
    if lam != 0:
        y_pred_cv = np.power((y_pred_boxcox * lam + 1), (1 / lam))
    else:
        y_pred_cv = np.exp(y_pred_boxcox)

    # Evaluate the model using regression metrics on the original scale
    mse = mean_squared_error(y, y_pred_cv)
    mae = mean_absolute_error(y, y_pred_cv)
    r2 = r2_score(y, y_pred_cv)

    # Print out the metrics
    print(f'Mean Squared Error (MSE): {mse:.4f}')
    print(f'Mean Absolute Error (MAE): {mae:.4f}')
    print(f'R² Score: {r2:.4f}')

    # Residual plot
    residuals = y - y_pred_cv
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=y_pred_cv, y=residuals)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Residual Plot')
    plt.grid(True)
    plt.show()

    # Predicted vs Actual plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=y, y=y_pred_cv)
    plt.plot([min(y), max(y)], [min(y), max(y)], color='red', linestyle='--')  # Line y=x
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')
    plt.title('Predicted vs Actual Values')
    plt.grid(True)
    plt.show()

    # Train the model on the full dataset for feature importance
    lasso_model.fit(X_scaled, y_scaled)

    # Get feature importance (coefficients from Lasso)
    importance = np.abs(lasso_model.coef_)

    # Create a DataFrame for the feature importance
    feature_importance_df = pd.DataFrame({
        'Feature': X_dummies.columns,
        'Importance': importance
    })

    # Sort by importance
    feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)

    # Plot feature importance
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Importance', y='Feature', data=feature_importance_df)
    plt.title('Feature Importance - Lasso Regression')
    plt.tight_layout()
    plt.show()