# %%
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score , mean_absolute_error, mean_squared_error , mean_absolute_percentage_error, median_absolute_error, mean_squared_log_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, BayesianRidge, HuberRegressor, PassiveAggressiveRegressor, RANSACRegressor, ARDRegression, OrthogonalMatchingPursuit
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.svm import SVR
from sklearn.linear_model import SGDRegressor
from sklearn.neural_network import MLPRegressor
import pickle


# %%
def get_data_separated(dataframe, test_size=0.2, random_state=None,choice=0):
    if choice == 0 :
        target_feature = ['ΔVC -p'] #+ ΔVC-m,ΔVC -p
    elif choice == 1 : 
        target_feature = ['ΔVC-m']
    else:
        target_feature = ['ΔVC-m','ΔVC -p']

    features_to_drop = ['substituent','molecular_formula','Name','canonical_smiles','cid']

    # Assuming 'target_feature' is the dependent variable (y)
    y = dataframe[target_feature]
    X = dataframe.drop(columns=features_to_drop + target_feature)
    X =X.fillna(0)

    # Split into training and testing sets, excluding data_val
    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )
    # create the next sets , you split the temp data into halfs 
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.5, random_state=random_state
    )

    return X_train, X_val, X_test, y_train, y_val, y_test

# %%
def scaling_data(X_train,X_val,X_test):

    scaler = StandardScaler()

    scaler.fit(X_train)
    
    # filename = '../models/scaler_for_all_features.sav'
    # pickle.dump(scaler, open(filename, 'wb'))
    
    scaled_train = scaler.transform(X_train)
    scaled_test = scaler.transform(X_test)
    scaled_val = scaler.transform(X_val)
    return scaled_train,scaled_val,scaled_test

# %%
def save_metrics_results(model,X_test,y_test,tag):

    y_pred = model.predict(X_test)
    r2 = r2_score(y_pred=y_pred,y_true=y_test)
    MAE = mean_absolute_error(y_pred=y_pred,y_true=y_test) 
    MSE = mean_squared_error(y_pred=y_pred,y_true=y_test) 
    MAPE = mean_absolute_percentage_error(y_pred=y_pred,y_true=y_test)
    MedAE = median_absolute_error(y_pred=y_pred,y_true=y_test)
    # MSLE = mean_squared_log_error(y_pred=y_pred,y_true=y_test)

    results = {
        "R^2": r2,
        "MAE":MAE,
        "MSE":MSE,
        "MAPE":MAPE ,
        "MedAE":MedAE,
        # "MSLE":MSLE
    }
    mean_df = pd.DataFrame(data=results,index=[f"{tag}"])
    return mean_df

def load_regression_models():
    models = {
        'Linear Regression': LinearRegression(),
        'Ridge': Ridge(),
        'Lasso': Lasso(),
        'Elastic Net': ElasticNet(),
        'Bayesian Ridge': BayesianRidge(),
        'Huber Regressor': HuberRegressor(),
        'Passive Aggressive Regressor': PassiveAggressiveRegressor(),
        'RANSAC Regressor': RANSACRegressor(),
        'ARD Regression': ARDRegression(),
        'Orthogonal Matching Pursuit': OrthogonalMatchingPursuit(),
        'Decision Tree': DecisionTreeRegressor(),
        'Random Forest': RandomForestRegressor(),
        'Gradient Boosting': GradientBoostingRegressor(),
        'K-Nearest Neighbors': KNeighborsRegressor(),
        'Support Vector Regression': SVR(),
        'Stochastic Gradient Descent': SGDRegressor(),
        'Multi-layer Perceptron': MLPRegressor()
    }
    return models



# %%
def plot_the_r2(y_true, X_pred , model,filename): #(model,X_test,y_test,tag)
    y_true = y_true.to_numpy().ravel()
    y_pred = model.predict(X_pred)
    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.scatter(y_true, y_pred, color='blue')
    plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], linestyle='--', color='red')
    plt.title('Real Values vs. Predicted Values')
    plt.xlabel('Actual Values')
    plt.ylabel('Predicted Values')

    # Plot Residuals
    residuals = y_true - y_pred

    plt.subplot(1, 2, 2)
    plt.scatter(y_true, residuals, color='green')
    plt.axhline(y=0, linestyle='--', color='red')
    plt.title('Residual Plot')
    plt.xlabel('Actual Values')
    plt.ylabel('Residuals')

    plt.tight_layout()
    plt.savefig(fname=filename ,transparent=True ,format='jpg')

    # plt.show()


# %% [markdown]
# * create a loop that goes over all the list of algoritms 
# * then use the fit function to train them
# * the next step is to save the metrics per model 
# * check the evaluation and plot the models that you think is relevant


# %%
def plot_features_using_permutation_importance(model,X_train, y_train,random_state, feature_names_ori,filename):
    feature_importances = permutation_importance( model, X_train, y_train, n_jobs=-1,random_state=random_state)
    sorted_idx = feature_importances.importances_mean.argsort()
    feature_names = [feature_names_ori[i] for i in sorted_idx]
    name = filename.split("/")[-1]
    plt.figure(figsize=(10, 6))
    plt.bar(feature_names, feature_importances.importances_mean[sorted_idx])
    plt.grid(True,which="major")
    plt.xlabel('Features')
    plt.ylabel('Importance')
    plt.title(f'Feature Importance {name}')
    plt.xticks(rotation=90)  # Rotate feature names for readability
    plt.tight_layout()
    # plt.show()
    plt.savefig(fname=filename ,transparent=True ,format='svg')
    return feature_names, feature_importances.importances_mean[sorted_idx]




# %% 
def store(b, file_name):
    pickle.dump(b, open(file_name, "wb"))

def load(file_name):
    b = {}
    try:
        b = pickle.load(open(file_name, "rb"))
        print("Loading Successful")
        return b
    except (OSError, IOError) as e:
        print("Loading Failed. Initializing to empty")
        b = {}
        return b
