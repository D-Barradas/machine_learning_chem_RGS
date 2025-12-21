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

Quick Start — SELFIES transfer-learning
--------------------------------------

If you want to reproduce the SELFIES-based transfer-learning pipeline (preprocessing → HF dataset → fine-tune), follow these quick steps.

1. Create the environment (recommended — includes RDKit):

    ```bash
    conda env create -f environment_selfies.yml
    conda activate selfies-chem-ml
    ```

    Or install pip dependencies (RDKit via conda recommended):

    ```bash
    python -m pip install --upgrade pip
    python -m pip install -r requirements.txt --user
    ```

2. Build the small HuggingFace dataset (uses propagated SELFIES):

    ```bash
    python src/build_hf_dataset.py --input data/substituent_selfies_propagated.csv \
         --output-prefix data/hf_dataset_small --mode small
    ```

    This writes `data/hf_dataset_small.csv` and the splits
    `data/hf_dataset_small_train.csv`, `_val.csv`, `_test.csv`.

3. (Optional) Expand the dataset from the large PubChem-derived table and convert SMILES→SELFIES:

    ```bash
    # make sure `selfies` is installed in the environment
    python src/build_hf_dataset.py --input data/dv_values_by_cid.csv \
         --output-prefix data/hf_dataset_large --mode large
    ```

4. Prepare training environment and run fine-tuning (example placeholder):

    - Install `transformers`, `datasets`, `accelerate` (see `requirements.txt`).
    - Use a SELFormer / SELFIES model checkpoint (or a compatible tokeniser) and fine-tune with a regression head.

5. Helpful utilities in `src/`:

    - `src/prepare_dataset_selfies.py` — convert SMILES → SELFIES and collect DV labels.
    - `src/map_names_to_smiles.py` — map substituent names to SMILES using PubChem.
    - `src/build_hf_dataset.py` — build HF-style CSVs and train/val/test splits.

Phenanthroline Derivative Generator
----------------------------------

Generate orthophenanthroline derivatives (positions 2–9) and score interaction energies using DV descriptors.

- Script: `src/generate_phen_derivatives.py`
- Scorer: `src/score_derivatives_eq8_like.py`
- DV mapping: `data/substituent_dv_mappings.csv` (aggregated from workspace data)

Usage examples

```bash
# Build DV mapping (once)
python src/build_substituent_library.py

# Generate mono/di derivatives (default substituents)
python src/generate_phen_derivatives.py \
    --out data/phen_derivatives_generated_fix.csv \
    --modes mono di

# Generate mono/di/tri/tetra using identical-only combos to avoid explosion
python src/generate_phen_derivatives.py \
    --out data/phen_derivatives_generated_fix.csv \
    --subs NO2 CN CHO OCF3 Cl Br NH2 OMe Me Et SiH3 \
    --modes mono di tri tetra \
    --identical-only

# Score derivatives using Eq.8-like linear fits (Table 3) and tri interpolation
python src/score_derivatives_eq8_like.py
```

Flags

- `--modes`: one or more of `mono`, `di`, `tri`, `tetra`.
- `--subs`: substituent names to attach; defaults to a core set.
- `--identical-only`: generate only identical-substituent combinations for given multiplicity.
- `--eint`: optional path to `data/e_int_consolidated.csv` to enrich DV_C mapping.

Scoring assumptions

- Per multiplicity (1, 2, 4), fit `E_int ≈ a + b · DV_C` using Table 3 (`Table3_eq8`).
- `DV_C` per derivative is the sum of substituent DV_C values (`dv_c_sum`).
- DV_C priority: `data/substituent_dv_mappings.csv` → Table 3 values → small defaults.
- Tri (3) predictions are interpolated between the fitted di (2) and tetra (4) models.

Outputs

- Generated: `data/phen_derivatives_generated_fix.csv` with `system`, `smiles`, `sites`, `substituents`, `n_subs`, `dv_c_sum`.
- Scored: `data/phen_derivatives_scored.csv` with `e_int_pred_eq8_fit` added.

Caveats

- Position labels (2–9) reflect RDKit aromatic C–H detection order (stable but approximate).
- Mixed tri/tetra combos can be enabled (remove `--identical-only`), but row counts grow combinatorially.

Phase 4: Reinforcement Learning for Ligand Optimization
------------------------------------------------------

Generate optimized phenanthroline derivatives using PPO over SELFIES tokens.

Components
- Reward service: `src/reward_service.py` wraps trained ExtraTrees predictor; reward = −E_int_pred + penalties.
- Tokenizer: `src/selfies_tokenizer.py` builds SELFIES vocab from SMILES data.
- RL trainer: `src/rl_optimize_phen_ppo.py` implements PPO with LSTM policy over SELFIES tokens.

Usage

```bash
# Train PPO policy to optimize predicted E_int
python src/rl_optimize_phen_ppo.py \
  --model models/eint_hpo_rf.pkl \
  --dv-mapping data/substituent_dv_mappings.csv \
  --smiles-data data/phen_derivatives_scored.csv \
  --out-dir results/rl_ppo \
  --episodes 200 \
  --batch-size 32 \
  --max-len 50 \
  --lr 1e-4
```

Outputs
- Trained policy: `results/rl_ppo/best_policy.pt`
- Top molecules: `results/rl_ppo/top_molecules.csv` with SMILES/SELFIES/reward/eint_pred
- Training logs: `results/rl_ppo/training_log.json` (mean/max reward per episode)

Notes
- Uses CPU by default (cuDNN compatibility issues on some systems).
- Reward combines predicted E_int (lower is better) with complexity penalties.
- Policy is LSTM over SELFIES tokens with KL regularization to base policy.
- Current surrogate uses Eq8 pseudo-labels; retrain on experimental E_int when available.

If you want, I can add a ready-to-run `train.py` (HuggingFace training loop) and a short notebook showing a complete PoC fine-tune on the small dataset.

--------


**PoC Outputs**

- **Script**: `src/selfies_tokenizer_poc.py` — builds a bracket-token SELFIES vocab, trains a tiny PyTorch regressor, and saves a checkpoint.
- **Primary model checkpoint**: `models/poc_selfies/poc_selfies.pth` — a PyTorch checkpoint dict containing `model_state_dict`, `stoi`, `itos`, and `args`.
- **Input files expected**: `data/hf_dataset_small_train.csv`, `data/hf_dataset_small_val.csv`, `data/hf_dataset_small_test.csv` (CSV with at least `text` and the target column such as `DV_C`).
- **Vocab format**: `stoi` (token -> id) and `itos` (id -> token) are stored inside the checkpoint; tokens are SELFIES bracket tokens (e.g., `[C]`, `[=O]`) plus special tokens `<PAD>` and `<UNK>`.
- **How to run (example)**:

```bash
python src/selfies_tokenizer_poc.py \
    --train data/hf_dataset_small_train.csv \
    --val data/hf_dataset_small_val.csv \
    --test data/hf_dataset_small_test.csv \
    --target DV_C \
    --out models/poc_selfies/poc_selfies.pth \
    --epochs 5 \
    --batch_size 2
```

- **Notes**: The PoC is intentionally small — use the larger HF dataset (build from `data/dv_values_by_cid.csv`) for meaningful training. The script falls back to a regex-based tokenizer if the `selfies` package is not installed.


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
- [Dr. Remya Nair](https://github.com/remya-gs)

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
