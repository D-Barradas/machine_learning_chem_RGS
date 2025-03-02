Regression CompChem 
==============================

## Computational Chemistry Project

This project focuses on querying the PubChem database to retrieve as many descriptors as possible for various chemical subtituents in aromatic corganic compounds. Additionally, it includes several engineered features to enhance the dataset. The project is organized into 12 Jupyter notebooks, each addressing different aspects of the analysis and machine learning processes:

1. **00_chem_comp_basic.ipynb**: Introduction to basic computational chemistry concepts.
2. **01_EDA_data.ipynb**: Exploratory Data Analysis (EDA) of the retrieved data.
3. **02_get_properties_per_name.ipynb**: Script to get properties based on chemical names.
4. **03_get_the_properties_per_line.ipynb**: Script to get properties based on line input.
5. **04_Vanilla_ML_approach.ipynb**: Initial machine learning approach to predict properties.
6. **05_ML_search_for_best_algorithm.ipynb**: Search for the best regression algorithm.
7. **06_HPO_for_best_model.ipynb**: Hyperparameter optimization for the best model obtained from NB 05.
8. **07_prototype_for_serial_HPO.ipynb**: Prototype for hyperparameter optimization for a series of regression algorithm.
9. **08_ML_search_best_algortim_multiple_cols.ipynb**: Search for the best algorithm for that predicts multiple values.
10. **09_serial_HPO_for_selected_algos_predic_delta_VC_m**: Serial HPO for selected algorithms to predict a single output.
11. **10_serial_HPO_for_selected_algos_predic_delta_VC_m_VC_p**: Serial HPO for selected algorithms to predict multiple outputs.
12. **11_Vanilla_ML_multipred_cols.ipynb**: Vanilla machine learning approach for prediction of a substituents multiple systems (Benzene, Pyrene, Pyridine, Butadiene, Butadiyne) columns .
13. **12_Vanilla_ML_multipred_multisubs.ipynb**: Vanilla machine learning approach for multiple substituents on multiple systems.

This comprehensive project leverages both standard and advanced machine learning techniques to analyze and predict properties from small substituents, providing a robust material for computational chemistry research.

--------


## Project Organization
------------

    ├── LICENSE
    │
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── data               <- A folder for the expected csv files for input and destinations for output 
    │    └── search_dataframes_result <- CSV files from the HPO for the models
    │
    ├── models             <- selected models after HPO 
    │
    ├── notebook          <- Jupyter notebooks from the projecs
    │
    ├── environment.yml   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `conda env export --no-builds | head -n -1 > environment.yml`
    │
    ├── figures            <- A folder for the expected figures for diverse results
    │                        
    ├── src                <- source code for use in this project.
    

--------

## Contributing

Contributions are welcome! Please follow these steps to contribute:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/your-feature-name`).
3. Commit your changes (`git commit -m 'Add some feature'`).
4. Push to the branch (`git push origin feature/your-feature-name`).
5. Open a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) for the chemical data.
- [cookiecutter data science project template](https://drivendata.github.io/cookiecutter-data-science/)

## Authors

- [Dr. Didier Barradas Bautista](https://www.github.com/D-barradas)
- [Dr. Remya Nair]())

## 🔗 Links

[KAUST Core Labs](https://corelabs.kaust.edu.sa/
) : 
[![linkedin](https://img.shields.io/badge/linkedin-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/company/kaust-core-labs/about/) [![twitter](https://img.shields.io/badge/twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/kaust_corelabs)


[KAUST Vizualization Core Lab](https://corelabs.kaust.edu.sa/labs/detail/visualization-core-lab) :
[![KVL](https://img.shields.io/badge/twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/KAUST_Vislab)  
[![YouTube Channel Views](https://img.shields.io/youtube/channel/views/UCR1RFwgvADo5CutK0LnZRrw?style=social)](https://www.youtube.com/channel/UCR1RFwgvADo5CutK0LnZRrw)


<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

<!-- datascience-kvl-template
==============================

Use this template to start you project on data science

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p> -->
