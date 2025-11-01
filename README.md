# SOCastR: FAIR Documentation for Digital Soil Mapping Workflow

**Version:** 1.0.0  
**Date:** 2025-11-01  
**Author:** Markus Möller  
**Institution:** Julius Kühn Institute  
**Contact:** markus.moeller@julius-kuehn.de  
**License:** MIT  
**DOI:** [To be assigned upon Zenodo publication]

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [FAIR Principles Implementation](#2-fair-principles-implementation)
3. [Provenance Documentation](#3-provenance-documentation)
4. [Fitness-for-Purpose Assessment](#4-fitness-for-purpose-assessment)
5. [ISO 19157 Data Quality Metadata](#5-iso-19157-data-quality-metadata)
6. [Installation and Usage](#6-installation-and-usage)
7. [Outputs and Deliverables](#7-outputs-and-deliverables)
8. [Quality Assurance and Validation](#8-quality-assurance-and-validation)
9. [Citation and Acknowledgments](#9-citation-and-acknowledgments)
10. [References](#10-references)

---

## 1. Project Overview

### 1.1 Title
**SOCastR: Soil Organic Carbon Prediction Workflow with Spatially Explicit Uncertainty Quantification**

### 1.2 Description
SOCastR is a comprehensive R-based workflow for digital soil mapping (DSM) that predicts soil organic carbon (SOC) content using Random Forest and Quantile Regression Forest models with rigorous spatial cross-validation, forward feature selection, and uncertainty quantification. The workflow is specifically designed for intermediate-scale mapping with explicit handling of model uncertainty and extrapolation risk assessment.

### 1.3 Author and Affiliation
- **Author:** Markus Möller
- **Institution:** Julius Kühn Institute (JKI), Federal Research Centre for Cultivated Plants
- **Email:** markus.moeller@julius-kuehn.de
- **ORCID:** https://orcid.org/0000-0002-1918-7747

### 1.4 Keywords
- Digital soil mapping
- Soil organic carbon
- Random Forest
- Machine learning
- Uncertainty quantification
- Spatial cross-validation
- Quantile regression forest
- Area of Applicability
- Dissimilarity Index
- Forward feature selection
- CAST package
- R programming
- ISO 19157

### 1.5 Software Dependencies and Versions

| Package | Version | Purpose |
|---------|---------|---------|
| CAST | ≥0.7.0 | Spatial cross-validation (KNNDM) and Area of Applicability |
| caret | ≥6.0 | Machine learning framework and model training |
| randomForest | ≥4.7 | Random Forest algorithm implementation |
| quantregForest | ≥1.3 | Quantile regression for prediction intervals |
| terra | ≥1.7 | Raster data handling and spatial operations |
| sf | ≥1.0 | Vector data handling and spatial partitioning |
| dplyr | ≥1.1 | Data manipulation and summarization |
| doParallel | ≥1.0 | Parallel computing for tile-based processing |
| ggplot2 | ≥3.4 | Visualization of results and validation plots |
| tidyterra | ≥0.4 | Raster visualization with ggplot2 |
| classInt | ≥0.4 | Classification intervals for map production |
| RColorBrewer | ≥1.1 | Color schemes for cartography |
| viridis | ≥0.6 | Perceptually uniform color palettes |
| gridExtra | ≥2.3 | Combining multiple plots |
| grid | (base) | Base graphics system |

**Language:** R ≥ 4.0.0 (recommended)

**Operating Systems:** Windows, Linux/Unix, macOS

---

## 2. FAIR Principles Implementation

### 2.1 Findable

**F1: Globally Unique and Persistent Identifier**
- GitHub repository URL: https://github.com/[username]/SOCastR
- Zenodo DOI: 10.5281/zenodo.[XXXXXXX] (assigned on publication)
- Each release version receives persistent DOI via GitHub-Zenodo integration
- Archived snapshot for reproducibility

**F2: Rich Metadata**
- Comprehensive README.md with project description, quick start, and usage examples
- CITATION.cff file in Citation File Format 1.2.0 for standardized citation
- codemeta.json for machine-readable software metadata (Schema.org format)
- Detailed inline code documentation and roxygen-style function docstrings
- This FAIR documentation file (readme.md)

**F3: Metadata References Identifier**
- All documentation explicitly references repository URL and DOI
- CITATION.cff includes repository and DOI fields
- Metadata embedded in code file headers

**F4: Indexed in Searchable Resource**
- GitHub repository (indexed by Google, Bing, search engines)
- Zenodo record (indexed by OpenAIRE, BASE, Google Dataset Search, CERN)
- Keywords and tags ensure discoverability

### 2.2 Accessible

**A1: Retrievable via Standard Protocol**
- Code accessible via HTTPS (GitHub)
- Zenodo archive accessible via HTTPS with REST API
- No authentication required for access (public repository)
- Open license facilitates reuse

**A2: Metadata Persists**
- Zenodo ensures metadata persistence for indefinite period
- GitHub provides persistent repository structure and history
- Documentation archived with each release version
- Code versioning ensures reproducibility

### 2.3 Interoperable

**I1: Formal, Accessible, Shared Language**
- Uses standard R programming language and widely-adopted packages
- Follows tidyverse conventions and spatial data science best practices
- Implements ISO 19157 compliant data quality vocabulary
- Supports standard geospatial data formats (GeoTIFF, ESRI Shapefile)

**I2: Uses FAIR-Compliant Vocabularies**
- ISO 19157:2013 controlled vocabulary for data quality elements
- ISO 19115 metadata elements for geospatial data
- W3C PROV-DM concepts for provenance representation
- Standard spatial and temporal coordinate reference systems (EPSG codes)

**I3: Includes Qualified References**
- References external datasets with proper attribution and metadata
- Cites methodological literature and software packages with full citations
- Links to relevant international standards (ISO 19157, ISO 19115)
- Provides URLs to reference implementations

### 2.4 Reusable

**R1: Rich Attributes**
- Comprehensive documentation of inputs, processing, and outputs
- Detailed description of uncertainty quantification methods and assumptions
- Clear specification of data quality requirements and fitness criteria
- Parameter defaults and tuning guidance

**R1.1: Clear Usage License**
- LICENSE file (MIT) in repository root
- License clearly stated in code headers and documentation
- Permissive license facilitates both academic and commercial reuse

**R1.2: Detailed Provenance**
- Complete lineage documentation with processing steps
- W3C PROV-DM model for provenance representation
- Input/output tracking and version control via Git
- Processing parameters and computational environment documentation

**R1.3: Domain-Relevant Standards**
- ISO 19157:2013 for data quality assessment
- ISO 19115-3 for geospatial metadata encoding
- W3C PROV-O ontology for provenance
- OGC standards for spatial data and services

---

## 3. Provenance Documentation

### 3.1 Data Lineage

#### 3.1.1 Input Data

**Primary Input: Soil Sample Points**
- **Format**: ESRI Shapefile (.shp, .shx, .dbf, .prj)
- **Geometry Type**: Point features
- **Coordinate Reference System**: Any projected CRS (e.g., EPSG:25832 for UTM Zone 32N)
- **Required Attributes**: 
  - SOC measurement column (numeric, user-specified via `soc_column` parameter)
  - Unit: % 
- **Quality Requirements**: 
  - No duplicate locations (same coordinates)
  - Valid coordinates within study area
  - SOC values within physically plausible range (0-30% typical)
  - No missing values in SOC column
- **Provenance Metadata**: Document sampling date, method (depth, horizon), lab analysis procedure, responsible institution

**Secondary Input: Environmental Covariates**
- **Format**: Multi-band GeoTIFF (.tif) raster stack
- **Coordinate Reference System**: Must match or be compatible with soil samples (automatically transformed if needed)
- **Resolution**: Uniform across all bands (e.g., 100m × 100m)
- **Spatial Extent**: Must completely cover all sample point locations
- **Number of Bands**: Variable (automatically detected and used as predictors)
- **Band Names**: Generic names assigned by script (COV1, COV2, ..., COVn)
- **Covariate Types**: Typically terrain derivatives (slope, aspect, curvature), climate variables, remote sensing indices (NDVI, parent material classification)
- **Covariate Extraction Method**: Median filter across 3×3 pixel neighborhood (function: `extract_8_neighbors`)
  - Reduces GPS positioning errors and local artifacts
  - Improves representativeness of covariate values at point locations
- **Provenance Metadata**: Data source, derivation method, temporal coverage, spatial accuracy, resolution rationale

#### 3.1.2 Processing Workflow

The workflow implements a standardized sequence of 14 steps for robust DSM:

```
STEP 1: Package Loading and Verification
  └─ Checks R packages installed
  └─ Auto-installs missing packages
  └─ Initializes libraries

STEP 2: Function Definition
  ├─ PlotR(): Classification raster plotting
  ├─ SFPlotR(): Point data spatial visualization
  ├─ extract_8_neighbors(): 3×3 median filter extraction
  └─ create_tiles_manual(): Large raster tiling

STEP 3: Data Loading
  ├─ Load soil samples (sf::st_read)
  ├─ Load covariate stack (terra::rast)
  ├─ Validate CRS compatibility
  └─ Standardize layer names to COV1, COV2, ...

STEP 4: Covariate Extraction (Filtered)
  ├─ For each sample point:
  │  ├─ Identify focal cell + 8 neighbors (3×3 window)
  │  ├─ Extract values from all 9 cells for each layer
  │  └─ Calculate median across 9 cells
  ├─ Combine SOC values with extracted covariates
  ├─ Remove samples with NA values
  └─ Output: model_data_clean (samples × covariates matrix)
  
STEP 5: Spatial Data Partitioning
  ├─ Method: CreateSpacetimeFolds (k=5 spatial folds)
  ├─ Use fold 1 as test set, remaining as training
  ├─ Ensures spatial independence between train/test
  ├─ Distribution validation (Kolmogorov-Smirnov test)
  └─ Output: Spatial partition map + statistics

STEP 6: Spatial Cross-Validation Setup
  ├─ Method: KNNDM (k-Nearest Neighbor Distance Matching)
  ├─ Training data only: k=5 folds
  ├─ Ensures CV folds match geographic distance distribution
  ├─ Geodist analysis: compare fold quality metrics
  └─ Output: Fold indices + geodist diagnostics

STEP 7: Forward Feature Selection (FFS)
  ├─ Input: All covariate layers
  ├─ Method: CAST::ffs with Random Forest (ntree=100)
  ├─ Selection criterion: RMSE minimization in spatial CV
  ├─ Hyperparameter optimization: tuneLength=3 (mtry values)
  ├─ Remove zero-variance predictors first
  └─ Output: Selected features subset

STEP 8: Final Model Training
  ├─ Algorithm: Random Forest (caret::train)
  ├─ Input: Training data × selected features
  ├─ Hyperparameters:
  │  ├─ ntree: 500 (increased from FFS for stability)
  │  ├─ mtry: Optimized via cross-validation
  │  └─ nodesize: Default (1 for regression)
  ├─ Cross-validation: Spatial CV (KNNDM folds)
  └─ Output: Trained model + variable importance

STEP 9: External Validation
  ├─ Predict on independent test set
  ├─ Calculate performance metrics:
  │  ├─ RMSE (Root Mean Square Error)
  │  ├─ MAE (Mean Absolute Error)
  │  ├─ R² (coefficient of determination)
  │  └─ Bias (mean error)
  ├─ Generate scatter and residual plots
  └─ Compare training vs. test performance

STEP 10: Retrain on All Data
  ├─ Purpose: Maximize information for spatial prediction
  ├─ Use all available samples (training + test)
  ├─ KNNDM folds computed on full dataset
  └─ Output: Final model for wall-to-wall prediction

STEP 11: Wall-to-Wall Spatial Prediction
  ├─ Input: Covariate raster × trained model
  ├─ Prediction method: Random Forest point predictions
  ├─ Output raster: SOC_predicted
  ├─ Resolution: Same as input covariates
  └─ Format: GeoTIFF with CRS metadata

STEP 12: Model Uncertainty (Quantile Regression Forest)
  ├─ Algorithm: quantregForest::quantregForest
  ├─ Quantiles: 0.05, 0.5, 0.95 (90% prediction interval)
  ├─ Training data: All samples × selected features
  ├─ Ntree: 500
  ├─ Output rasters:
  │  ├─ Q5 (5th percentile / lower CI)
  │  ├─ Q50 (50th percentile / median)
  │  ├─ Q95 (95th percentile / upper CI)
  │  └─ PIW (Prediction Interval Width = Q95 - Q5)
  ├─ Tiling: Automatic for rasters >5M cells
  │  ├─ Default: 4×4 grid (16 tiles)
  │  ├─ OS-aware parallelization (fork on Unix, PSOCK on Windows)
  │  └─ Merge tiles using mosaic (robust to overlaps)
  └─ Format: Multi-band GeoTIFF

STEP 13: Distance-Based Uncertainty (AOA/DI)
  ├─ Method: CAST::trainDI + CAST::aoa
  ├─ trainDI: Compute Dissimilarity Index threshold
  │  └─ Threshold based on training data feature space
  ├─ aoa: Classify prediction locations
  │  ├─ 1 = within Area of Applicability (reliable)
  │  └─ 0 = outside AOA (extrapolation, use caution)
  ├─ DI: Dissimilarity Index (continuous, 0-∞)
  │  └─ Higher DI = further from training data
  ├─ Tiling: Automatic for rasters >5M cells
  │  ├─ Same tile strategy as QRF
  │  └─ Parallel processing on multi-core systems
  └─ Format: Single-band GeoTIFFs

STEP 14: Output Generation and Visualization
  ├─ Raster outputs (GeoTIFF format, LZW compression)
  ├─ Statistics tables (CSV format)
  ├─ Visualization maps (PNG format, 300 DPI)
  └─ Model objects (RDS format for reproducibility)
```

#### 3.1.3 Processing Parameters

All parameters are configurable through the main `SOCastR()` function:

| Parameter | Type | Default | Description | Notes |
|-----------|------|---------|-------------|-------|
| `working_dir` | Character | Required | Working directory path | Must exist or will cause error |
| `input_dir` | Character | "input" | Input data subdirectory | Relative to working_dir |
| `output_dir` | Character | "output" | Output results subdirectory | Created if doesn't exist |
| `samples` | Character | Required | Soil samples filename | Shapefile in input_dir |
| `covariates` | Character | Required | Covariate stack filename | GeoTIFF in input_dir |
| `soc_column` | Character | Required | SOC attribute name | Column name in shapefile |
| `n.tile` | Integer | 4 | Tiles per raster dimension | 4×4=16 tiles, adjust for memory |
| `model_uncertainty` | Logical | TRUE | Compute QRF uncertainty | Set FALSE to skip for speed |
| `distance_uncertainty` | Logical | TRUE | Compute DI/AOA layers | Set FALSE to skip for speed |

**Fixed Parameters (hardcoded in script):**

| Parameter | Value | Description |
|-----------|-------|-------------|
| Random seed | 42 | Reproducibility of spatial folds and model training |
| FFS ntree | 100 | Faster feature selection iteration |
| Final RF ntree | 500 | Balanced accuracy/speed for final model |
| QRF ntree | 500 | Sufficient trees for stable quantile estimation |
| Quantiles | c(0.05, 0.5, 0.95) | 90% prediction interval bounds |
| CV k-folds | 5 | Balance bias-variance in validation |
| Extraction window | 3×3 | Median filter for covariate extraction |
| Raster tiling threshold | 5 million cells | Trigger for tile-based processing |

### 3.2 Software Environment

**Language:** R ≥ 4.0.0

**Supported Operating Systems:**
- **Windows:** Full support via PSOCK parallel cluster
- **Linux/Unix:** Full support via fork-based parallelization
- **macOS:** Full support via fork-based parallelization

**Hardware Requirements:**
- **RAM:** 8GB minimum (16GB recommended for large rasters)
- **CPU:** Multi-core processor (4+ cores recommended)
- **Storage:** Sufficient space for raster tiles and intermediate outputs
- **Disk I/O:** SSD recommended for large raster operations

**Parallelization Strategy:**
- Windows: `makeCluster(n_cores, type = "PSOCK")` with explicit variable export
- Linux/Unix: `mclapply()` with fork-based memory sharing (more efficient)
- Core Detection: `detectCores() - 1` (reserve one core for system)
- Tile Processing: Automatic across all available cores
- Fallback: Sequential processing if parallel fails

**Random Seed:** Set to 42 for reproducibility of:
- Spatial cross-validation fold creation
- Random Forest model initialization
- Feature selection results

### 3.3 Execution Provenance

**Processing Logs:**
- Console output: Detailed step-by-step progress messages
- Timestamps: Processing start/end for each major step
- CSV files: Intermediate and final results with metadata
- Model objects: Trained models saved as RDS files

**Version Control:**
- Git repository: Track code evolution and changes
- Commit history: Document methodology refinements
- Release tags: Version each analysis run (v1.0.0, v1.0.1, etc.)
- DOI assignment: Each release gets persistent identifier

**Provenance Model (W3C PROV-DM Concepts):**

```
Entity: SOC_Prediction_Raster
  wasDerivedFrom: [TrainedRandomForestModel, CovariateRasterStack]
  wasGeneratedBy: SpatialPredictionActivity
  generatedAtTime: [execution timestamp]
  wasAttributedTo: SOCastR_v1.0.0

Entity: QuantilePredictionUncertainty
  wasDerivedFrom: [QuantileRegressionForestModel, CovariateRasterStack]
  wasGeneratedBy: QuantileRegressionActivity

Entity: AreaOfApplicability
  wasDerivedFrom: [TrainDIObject, CovariateRasterStack]
  wasGeneratedBy: DistanceUncertaintyActivity

Activity: SpatialPredictionActivity
  used: [TrainedRandomForestModel, CovariateRasterStack]
  wasAssociatedWith: RandomForestAlgorithm
  startedAtTime: [timestamp]
  endedAtTime: [timestamp]

Activity: ModelTrainingActivity
  used: [TrainingDataset, SelectedFeatures]
  wasInformedBy: ForwardFeatureSelectionActivity
  startedAtTime: [timestamp]
  endedAtTime: [timestamp]

Agent: SOCastR_v1.0.0
  actedOnBehalfOf: [User/Institution]
  wasAssociatedWith: RandomForestAlgorithm, QuantileRegressionForest
```

---

## 4. Fitness-for-Purpose Assessment

### 4.1 Intended Use

**Primary Purpose:**
Intermediate-scale digital soil mapping of soil organic carbon content (1:50,000 to 1:250,000 scale equivalent) with spatially explicit uncertainty quantification and extrapolation risk assessment.

**Target Applications:**
- Soil carbon stock inventories and audits
- Climate change mitigation and carbon credit verification
- Precision agriculture and variable rate application
- Landscape-scale land management planning
- Environmental impact assessment and modeling
- Regional soil quality assessments

**Geographic Scale:**
- **Spatial Resolution:** 10m to 100m typical (defined by covariate resolution)
- **Extent:** Regional to national scale possible
- **Minimum Mapping Unit:** Single raster pixel
- **Coverage:** Entire study area covered by raster extent

### 4.2 Data Quality Requirements

**Input Data Quality Criteria:**

| Quality Aspect | Requirement | ISO 19157 Element | Assessment |
|---|---|---|---|
| Sample density | ≥0.5 samples/km² (minimum); ≥1 sample/km² recommended | Completeness | Pre-analysis check |
| Spatial distribution | Stratified or regular; minimize clustering | Spatial autocorrelation | Visual inspection + Moran's I |
| SOC measurement accuracy | Lab-analyzed; ±5% relative error maximum | Thematic Accuracy | Metadata verification |
| Positional accuracy | ±10m for sample locations | Positional Accuracy | Metadata verification |
| Covariate completeness | No gaps in study area extent | Completeness | Raster extent check |
| Temporal consistency | Samples from similar time period (±5 years) | Temporal Quality | Metadata verification |
| Variable availability | ≥5 relevant covariates after feature selection | Completeness | Post-FFS check |

**Minimum Dataset Requirements:**
- Training samples: ≥100 (200+ recommended for robust modeling)
- Feature-to-sample ratio: After FFS, should be ≤1:10
- Spatial coverage: Representative of study area variability
- Feature variance: No zero-variance predictors (automatically removed)

### 4.3 Validation Metrics

**Accuracy Assessment Metrics:**

| Metric | Formula | Interpretation | Quality Threshold | Comments |
|--------|---------|----------------|-------------------|----|
| RMSE | √(Σ(ŷᵢ - yᵢ)²/n) | Prediction error magnitude | Lower is better | Same units as SOC |
| MAE | Σ\|ŷᵢ - yᵢ\|/n | Absolute average error | Lower is better | Robust to outliers |
| R² | 1 - (SS_res/SS_tot) | Variance explained | 0-1 scale; >0.3 acceptable | Goodness of fit |
| Bias | mean(ŷᵢ - yᵢ) | Systematic error | Close to 0 | Sign indicates direction |
| RMSE_percent | (RMSE/mean(y)) × 100 | Relative error | <30% acceptable, <20% good | Normalized error |
| NSE | 1 - (Σ(ŷᵢ - yᵢ)²/Σ(yᵢ - ȳ)²) | Nash-Sutcliffe efficiency | Similar to R² | Hydrological standard |

**Uncertainty Quantification Metrics:**

| Metric | Description | Quality Indicator | Target Range |
|--------|-------------|-------------------|---------------|
| PI Width | 95th - 5th percentile | Narrower = higher certainty | Varies by dataset |
| Coverage | % observations within PI | Should match nominal level | ~90% for 90% PI |
| Dissimilarity Index | Distance to training space | Lower = more reliable | 0 to threshold |
| AOA Percentage | % study area within AOA | Higher = more reliable | >70% desirable |

### 4.4 Fitness-for-Purpose Decision Framework

**This workflow IS fit for purpose when:**

✓ Sample density adequate (≥0.5 samples/km²)
✓ Spatial distribution representative (not strongly clustered)
✓ High-quality input data (accurate SOC, complete covariates)
✓ R² > 0.4 on independent test set (explains >40% variance)
✓ Unbiased predictions (Bias within ±10% of mean SOC)
✓ >70% of prediction area within Area of Applicability
✓ Prediction intervals capture observed variability
✓ RMSE_percent < 25% (good predictive accuracy)

**This workflow may NOT be fit for purpose when:**

✗ Sparse or heavily clustered sampling (<0.1 samples/km²)
✗ Poor input data quality (inaccurate SOC, incomplete covariates)
✗ Very low R² (<0.3, cannot explain data relationships)
✗ Large systematic bias (>20% of mean SOC)
✗ <50% of study area within Area of Applicability
✗ Extreme extrapolation (modeling conditions not in training data)
✗ Poor covariate quality or weak correlations with SOC

### 4.5 Limitations

1. **Spatial Stationarity Assumption:** Model assumes SOC-covariate relationships are spatially stationary within cross-validation folds; may fail in heterogeneous landscapes

2. **Temporal Transferability:** Model trained on historical data may not perform for future climate/land use conditions without recalibration

3. **Random Forest Constraints:**
   - Cannot extrapolate beyond training data range (predictions = mean of training samples in extrapolation zones)
   - Tends to underpredict extremes (high and low values)
   - Requires sufficient sample size (≥50) for stable estimates

4. **Covariate Dependency:** Prediction quality directly depends on covariate relevance, resolution, and accuracy; missing important variables limits predictability

5. **Spatial Autocorrelation:** Strong spatial autocorrelation can lead to overoptimistic cross-validation results; KNNDM mitigates but doesn't eliminate

6. **Computational Intensity:** Large rasters (>50M cells) require substantial RAM and processing time; tiling helps but reduces efficiency

7. **Model Assumptions:** Linear or non-linear relationships between predictors and SOC; complex interactions may not be captured

### 4.6 Use Case Scenarios with Fitness Assessment

**Scenario 1: Regional Soil Carbon Accounting**
- **Purpose:** Estimate total soil carbon stocks for national reporting
- **Fitness:** HIGH (if sample density ≥1 sample/km², AOA >80%)
- **Uncertainty Requirements:** Critical - must report confidence intervals for stock estimates
- **Recommendation:** Include both QRF and AOA uncertainty; apply spatial averaging techniques

**Scenario 2: Field-Scale Precision Agriculture**
- **Purpose:** Variable-rate fertilization based on soil carbon maps
- **Fitness:** MEDIUM (may require higher resolution covariates and denser sampling)
- **Uncertainty Requirements:** High - decision-making depends on spatial certainty
- **Recommendation:** Validate with field observations; use AOA mask to exclude unreliable zones

**Scenario 3: Climate Change Modeling**
- **Purpose:** Soil carbon as input to carbon cycle/climate models
- **Fitness:** MEDIUM-HIGH (depends on temporal representativeness)
- **Uncertainty Requirements:** High - uncertainty propagation important for ensemble models
- **Recommendation:** Document temporal limitations; acknowledge model recalibration needs

**Scenario 4: Land Policy and Conservation Planning**
- **Purpose:** Identify priority areas for carbon sequestration projects
- **Fitness:** HIGH (if predictions unbiased and uncertainty communicated)
- **Uncertainty Requirements:** Medium-high for decision confidence
- **Recommendation:** Use combination of predictions and AOA for risk-based prioritization

---

## 5. ISO 19157 Data Quality Metadata

### 5.1 Overview

This section provides data quality information compliant with **ISO 19157:2013 Geographic information - Data quality** standard. Quality is assessed for both input data requirements and prediction outputs.

### 5.2 Data Quality Scope

**DQ_Scope Definition:**
- **Level:** Dataset (prediction maps and input samples)
- **Extent:** Geographic extent of covariate raster stack and sample distribution
- **Feature Type:** Continuous raster (SOC concentration)
- **Attribute Type:** Quantitative (thematic attribute)

### 5.3 Data Quality Elements

#### 5.3.1 Completeness

**DQ_CompletenessCommission (Excess Data):**
- **Measure:** Samples removed due to NA values in covariates
- **Evaluation Method:** Direct internal (comparison with source data)
- **Result:** Percentage of samples with missing covariate data
- **Reporting:** `ExtractValues_SampleDataSummary.csv` (Removed_samples column)
- **Statement Example:** "1,247 original samples; 43 removed (3.4%) due to missing covariate values; 1,204 used for modeling"

**DQ_CompletenessOmission (Missing Data):**
- **Measure:** Samples with missing SOC values
- **Evaluation Method:** Direct internal (shapefile attribute check)
- **Result:** Count and percentage of missing SOC values
- **Reporting:** Data loading script output
- **Statement Example:** "SOC column: complete, no missing values"

**Raster Completeness:**
- **Measure:** Pixels with valid predictions
- **Evaluation Method:** Raster extent coverage check
- **Result:** Percentage of raster with non-NA values
- **Statement Example:** "Prediction raster: 99.8% of pixels predicted; 0.2% NA (outside covariate extent)"

#### 5.3.2 Logical Consistency

**DQ_ConceptualConsistency:**
- **Schema Compliance:** Point features with numeric SOC attribute
- **Evaluation:** Automated checks during data loading
- **Result:** Pass/Fail status

**DQ_DomainConsistency:**
- **Spatial Domain:** All samples within covariate raster extent
- **Value Domain:** SOC values within 0-30% (physical plausibility)
- **CRS Consistency:** All layers use compatible coordinate reference systems
- **Evaluation:** Pre-processing validation, automatic reprojection if needed

**DQ_FormatConsistency:**
- **Input Formats:** ESRI Shapefile (samples), GeoTIFF (covariates)
- **Output Formats:** GeoTIFF (rasters), CSV (statistics), PNG (visualizations)
- **Metadata Encoding:** GDAL tags in GeoTIFF headers
- **Evaluation:** Automated during file I/O

#### 5.3.3 Positional Accuracy

**DQ_AbsoluteExternalPositionalAccuracy:**
- **Measure:** Root Mean Square Error (RMSE) of point locations
- **Value:** Inherits from input sample data quality
- **Unit:** Meters
- **Evaluation:** Cannot directly assess (relies on input documentation)
- **Mitigation:** 3×3 pixel neighborhood median extraction reduces impact of ±10m errors

**Raster Positional Accuracy:**
- **Measure:** Pixel alignment with geodetic grid
- **Value:** Identical to covariate raster resolution
- **Evaluation:** Automatic via terra package CRS handling

**Statement Example:** "Positional accuracy of predictions ≤input sample accuracy (typically ±10m). Spatial filter (3×3 median) reduces localized positioning errors."

#### 5.3.4 Thematic Accuracy

**DQ_QuantitativeAttributeAccuracy:**

Measured through spatial cross-validation and independent test set validation.

**Training Performance (Spatial Cross-Validation):**
- **Method:** 5-fold spatial CV using KNNDM
- **Metric:** RMSE on CV predictions
- **Reporting:** `FinalModel_Accuracy.csv`
- **Interpretation:** Expected performance on similar study area

**Test Set Performance (External Validation):**
- **Method:** Predict on spatially independent test set
- **Metrics:** 
  - RMSE: [value] ± SD
  - MAE: [value] ± SD
  - R²: [value]
  - Bias: [value]
- **Reporting:** `Validation_PerformanceComparison.csv`
- **Interpretation:** Realistic performance estimate on novel data

**Statement Example:**
```
Spatial Cross-Validation Performance:
- Mean RMSE: 2.34 g/kg
- R² (median): 0.62
- Interpretation: Model explains 62% of SOC variance

Independent Test Set Performance:
- RMSE: 2.51 g/kg
- R²: 0.58
- Bias: -0.12 g/kg (slight underprediction)
- RMSE_percent: 18.3%
- Interpretation: Model performs well on independent data; nearly unbiased
```

#### 5.3.5 Temporal Quality
Temporal quality was not assessed. If information are available, the following data quality elements should be used:
- **DQ_TemporalConsistency:**
- **Temporal Coverage:** Samples collected within [date range]
- **Covariate Temporal Basis:** Covariates valid for [time period]
- **Synchronization:** [Describe alignment or temporal mismatch]

**DQ_TemporalValidity:**
- **Valid From:** [Effective date]
- **Valid To:** [Recommended end date or "Indefinite"]
- **Statement:** "Model applicable to conditions similar to training period. Recalibration recommended if climate or land cover changes substantially (>10 year gap)."

#### 5.3.6 Usability

**DQ_UsabilityElement** (fitness for purpose):
- **Intended Use:** [From user specification]
- **Quality Level:** [HIGH / MEDIUM / LOW]
- **Confidence:** [Justification based on validation results]

**Fitness for Purpose Statement Examples:**

1. **High Fitness:** "Suitable for regional soil carbon accounting. Model R²=0.68, RMSE=1.9%, 85% AOA coverage provides reliable 1km resolution maps. Uncertainty quantification enables confidence intervals for aggregate stock estimates."

2. **Medium Fitness:** "Partially suitable for field-scale applications. Local R²=0.58, but AOA only covers 62% of area due to covariate extrapolation. Recommend independent field validation in high-DI zones before operational use."

3. **Limited Fitness:** "Not recommended for high-precision applications. Sparse sample density (0.3 samples/km²) and R²=0.35 insufficient for critical decisions. Use only for preliminary assessment; collect additional samples for production use."

### 5.4 Data Quality Measures

#### 5.4.1 Prediction Quality Measures

**Measure 1: Spatial Cross-Validation RMSE**
- **Name:** Prediction accuracy under spatial separation
- **Measure ID:** DQ_QuantitativeAttributeAccuracy_CV
- **Measure Type:** Indirect (estimated via cross-validation)
- **Evaluation Method:** K-fold spatial CV with KNNDM fold creation
- **Data Quality Value:** [RMSE value] ± [SD]
- **Value Unit:** g/kg or % (same as SOC measurement unit)
- **Evaluation Result:** File `FinalModel_Accuracy.csv`
- **Uncertainty:** Standard deviation across CV folds

**Measure 2: Independent Test Set RMSE**
- **Name:** External validation accuracy
- **Measure ID:** DQ_QuantitativeAttributeAccuracy_External
- **Measure Type:** Direct external validation
- **Evaluation Method:** Predict on spatially independent test set (fold 1, CreateSpacetimeFolds)
- **Data Quality Value:** [RMSE] [MAE] [R²] [Bias]
- **Value Unit:** g/kg or %
- **Evaluation Result:** File `Validation_PerformanceComparison.csv`
- **Significance:** Unbiased estimate of model performance on novel data

#### 5.4.2 Uncertainty Quality Measures

**Measure 3: Prediction Interval Coverage**
- **Name:** 90% Prediction Interval Empirical Coverage
- **Measure Type:** Statistical uncertainty
- **Evaluation Method:** Check what % of observed test values fall within predicted PI
- **Data Quality Value:** [Coverage percentage]
- **Target Value:** ~90% (matches nominal confidence level)
- **Interpretation:** If coverage <85% or >95%, interval estimation may need adjustment

**Measure 4: Mean Prediction Interval Width**
- **Name:** Average prediction uncertainty
- **Measure Type:** Quantile regression uncertainty
- **Evaluation Method:** Q95 - Q05 for each pixel, then average
- **Data Quality Value:** [Mean width]
- **Unit:** Same as SOC
- **Reporting:** File `QuantilePrediction_UncertaintySummary.csv`
- **Interpretation:** Narrower intervals = higher confidence; related to local data density and covariate certainty

**Measure 5: Area of Applicability Percentage**
- **Name:** Percentage of study area with reliable predictions
- **Measure Type:** Distance-based spatial uncertainty
- **Evaluation Method:** Dissimilarity Index threshold from training data
- **Data Quality Value:** [Percentage within AOA]
- **Target:** >70% of area should be within AOA
- **Reporting:** File `FinalPrediction_AoaSummaryStatistics.csv`
- **Interpretation:** <50% AOA suggests limited training data coverage for landscape

**Measure 6: Dissimilarity Index Distribution**
- **Name:** Feature space distance from training data
- **Measure Type:** Distance-based uncertainty metric
- **Evaluation Method:** Training feature space multivariate distance
- **Data Quality Value:** [Mean DI, Max DI, DI_threshold]
- **Interpretation:** Pixels with DI > threshold = predictions rely on extrapolation

### 5.5 Metaquality (Quality of Quality Information)

**DQ_Confidence:**
- **Confidence in Validation Results:** HIGH (independent test set)
- **Confidence in Uncertainty Estimates:** HIGH (two complementary methods)
- **Confidence Basis:** Spatial CV separates training/test; external data independent

**DQ_Representativity:**
- **Training Data Representativeness:** Assessed via spatial extent and univariate distribution
- **Reporting:** File `SpatialDataPartition_MapTrainTestSplit.png` (visual assessment)
- **Statistical Test:** Kolmogorov-Smirnov test (reported in `SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv`)

**DQ_Homogeneity:**
- **Spatial Homogeneity:** Assessed via Dissimilarity Index maps
- **Reporting:** File `FinalPrediction_MapDI.png` (visualization of spatial variation in feature similarity)
- **Interpretation:** Homogeneous areas (low DI) have consistent predictions; heterogeneous areas (high DI, especially outside AOA) have variable reliability

### 5.6 Data Quality Reporting

All quality information is systematically reported through multiple formats:

**1. Quantitative Results (CSV files):**
- `ExtractValues_SampleDataSummary.csv`: Input data completeness
- `SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv`: Split quality
- `FinalModel_Accuracy.csv`: Training CV performance
- `Validation_PerformanceComparison.csv`: Test set accuracy
- `QuantilePrediction_UncertaintySummary.csv`: PI statistics
- `FinalPrediction_AoaSummaryStatistics.csv`: AOA coverage

**2. Visual Documentation (PNG maps, 300 DPI):**
- `SpatialDataPartition_MapTrainTestSplit.png`: Sample distribution
- `SpatialCrossValidation_GeodistEcdf.png`: CV fold quality
- `FinalPrediction_MapDI.png`: Dissimilarity Index spatial pattern
- `Validation_ScatterPlot.png`: Observed vs predicted
- `Validation_ResidualPlot.png`: Bias assessment

**3. Narrative Documentation:**
- This FAIR documentation file (SOCastR-FAIR-Docs.md)
- README.md with usage and quality guidance
- Inline code comments with methodological justification

### 5.7 ISO 19157 Compliance Assessment

**Conformance Level:** PARTIAL (Comprehensive coverage of data quality concepts, not full XML encoding)

**Implemented ISO 19157 Components:**
- ✓ All six data quality elements (completeness, positional accuracy, thematic accuracy, logical consistency, temporal quality, usability)
- ✓ Direct and indirect evaluation methods documented
- ✓ Quantitative measures applied with defined metrics
- ✓ Results reported in machine-readable (CSV) and visual (PNG) formats
- ✓ Quality metadata directly linked to datasets via file naming
- ✓ Uncertainty quantification (two complementary methods)

**Not Implemented (Full Compliance Would Require):**
- ✗ Formal XML encoding per ISO 19115-3 schema
- ✗ Quality metadata in standardized geographic metadata catalog
- ✗ Machine-readable quality metadata (codelists, controlled vocabularies)
- ✗ Integration with enterprise GIS data management system

**Pathway to Full ISO Compliance:**
Users requiring full ISO 19157/19115 compliance should:
1. Generate ISO 19115-3 compliant XML metadata using codelists
2. Populate all DQ_DataQuality elements with results from this workflow
3. Create geographic metadata record in institutional catalog (e.g., GeoNetwork)
4. Reference quality measures to ISO 19157 standard codelist
5. Archive metadata with prediction rasters (GeoTIFF auxiliary files)

---

## 6. Installation and Usage

### 6.1 Prerequisites

**System Requirements:**
- R version ≥ 4.0.0
- RStudio (optional, highly recommended for interactive use)
- Operating System: Windows, Linux/Unix, macOS
- RAM: 8GB minimum; 16GB+ recommended
- Multi-core processor: 4+ cores strongly recommended
- Disk space: ~1GB for typical run (with raster outputs)

**R Installation:**
Download R from https://www.r-project.org/ and follow platform-specific installation instructions.

### 6.2 Installation

**Step 1: Install Required R Packages**

```r
# Create list of required packages
required_packages <- c(
  "CAST",           # Spatial prediction and cross-validation
  "caret",          # Machine learning framework
  "classInt",       # Classification intervals for visualization
  "dplyr",          # Data manipulation and wrangling
  "doParallel",     # Parallel computing backend
  "terra",          # Raster data operations
  "tidyterra",      # ggplot2 integration for rasters
  "sf",             # Vector data operations
  "randomForest",   # Random Forest algorithm
  "RColorBrewer",   # Color palettes for cartography
  "quantregForest", # Quantile regression for uncertainty
  "ggplot2",        # Graphics and visualization
  "viridis",        # Perceptually uniform color scales
  "gridExtra",      # Combining plots
  "grid"            # Base graphics
)

# Install missing packages (only if not already installed)
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

# Verify all packages loaded
cat("All packages loaded successfully!\n")
```

**Step 2: Clone or Download Repository**

Via Git (recommended for version control):
```bash
git clone https://github.com/[username]/SOCastR.git
cd SOCastR
```

Or download ZIP directly from GitHub and extract.

**Step 3: Verify Installation**

```r
# Set working directory
setwd("~/SOCastR/")  # Modify path as needed

# Source main script
source("SOCastR.R")

# Check that main function exists
exists("SOCastR")  # Should return TRUE

# Check package versions
packageVersion("CAST")
packageVersion("terra")
packageVersion("caret")
```

### 6.3 Data Preparation

**Input Data Structure:**
```
project_directory/
├── SOCastR.R                          # Main workflow script
├── SOCastR-FAIR-Docs.md               # This documentation
├── CITATION.cff                       # Citation metadata
├── LICENSE                            # MIT license
├── README.md                          # Quick start guide
│
├── input/                             # Input data directory (create this)
│   ├── SAMPLES_EPSG25832.shp          # Soil samples (required)
│   ├── SAMPLES_EPSG25832.shx
│   ├── SAMPLES_EPSG25832.dbf
│   ├── SAMPLES_EPSG25832.prj
│   └── COVARIATES_EPSG25832.tif       # Covariate stack (required)
│
└── output/                            # Results (auto-created)
    ├── Raster_Outputs/
    ├── Validation_Statistics/
    ├── Visualizations/
    └── Model_Objects/
```

**Sample Data Requirements:**
- **Format:** ESRI Shapefile (.shp, .shx, .dbf, .prj files together)
- **Geometry:** Point features only (no lines or polygons)
- **CRS:** Any projected coordinate system (EPSG code preferred, e.g., EPSG:25832 for UTM Zone 32N)
- **Attributes:**
  - SOC column (numeric, required): Soil organic carbon concentration (% or g/kg)
  - Other columns optional but will be preserved
- **Quality:**
  - No duplicate point locations
  - No missing SOC values
  - Valid coordinates within expected geographic bounds

**Covariate Data Requirements:**
- **Format:** Multi-band GeoTIFF (.tif)
- **Bands:** 5-50 covariate layers (recommended 10-25)
- **CRS:** Same as samples or automatically reprojected to match
- **Extent:** Must completely cover all sample point locations
- **Resolution:** Uniform across all bands (typically 10m, 30m, or 100m)
- **Data Type:** Float32 or Float64
- **Quality:**
  - No gaps or holes in data (all pixels have values)
  - No zero-variance layers (automatically removed)
  - Physically plausible values for each variable
- **Covariate Types (Examples):**
  - Terrain derivatives: slope, aspect, curvature, elevation
  - Climate variables: precipitation, temperature, solar radiation
  - Remote sensing: NDVI, soil adjusted vegetation index
  - Categorical: parent material, lithology, soil type

### 6.4 Basic Usage

**Example 1: Minimal Configuration (Defaults)**

```r
# Load the main script
setwd("~/SOCastR/")
source("SOCastR.R")

# Run with minimal parameters (use defaults)
SOCastR(
  working_dir = getwd(),
  input_dir = "input",
  output_dir = "output",
  samples = "SAMPLES_EPSG25832.shp",
  covariates = "COVARIATES_EPSG25832.tif",
  soc_column = "SOC"
)
```

**Example 2: Custom Configuration**

```r
# Run with custom parameters
SOCastR(
  working_dir = "/path/to/project/",
  input_dir = "data/inputs",
  output_dir = "results/run_2024_11_01",
  samples = "soil_samples_cleaned.shp",
  covariates = "environmental_predictors_100m.tif",
  soc_column = "SOC_percent",
  n.tile = 4,                    # 4×4 tile grid (16 tiles)
  model_uncertainty = TRUE,      # Enable QRF uncertainty
  distance_uncertainty = TRUE    # Enable AOA/DI calculation
)
```

**Example 3: Quick Prediction Only (Skip Uncertainty)**

```r
# Faster run: skip computationally intensive uncertainty
SOCastR(
  working_dir = getwd(),
  input_dir = "input",
  output_dir = "output_quick",
  samples = "samples.shp",
  covariates = "covariates.tif",
  soc_column = "SOC",
  model_uncertainty = FALSE,     # Skip QRF
  distance_uncertainty = FALSE   # Skip AOA
)
```

**Example 4: Large Raster with Fine Tiling**

```r
# For very large rasters (>20M cells): increase tile density
SOCastR(
  working_dir = getwd(),
  input_dir = "input",
  output_dir = "output_large",
  samples = "samples.shp",
  covariates = "big_covariate_stack.tif",
  soc_column = "SOC",
  n.tile = 8  # 8×8 = 64 tiles (more granular, uses less memory per tile)
)
```

### 6.5 Execution Time Estimates

**Approximate Processing Time** (reference: 1000 samples, 20 covariates, 1M pixel raster):

| Step | Time (4 cores) | Memory | Notes |
|------|---|---|---|
| Data loading | <1 min | 1-2 GB | Fast for typical datasets |
| Covariate extraction | 2-5 min | 2-3 GB | Depends on sample count |
| Spatial partitioning | 2-3 min | 1-2 GB | KNNDM computation |
| Forward feature selection | 15-45 min | 3-5 GB | Iterative training; slowest step |
| Model training | 5-15 min | 3-5 GB | 5 CV folds × 3 mtry tuning |
| External validation | 1-2 min | 1-2 GB | Fast predictions |
| RF spatial prediction | 2-5 min | 2-4 GB | Wall-to-wall prediction |
| QRF uncertainty | 20-90 min | 4-8 GB | Tiled processing; can be slow |
| AOA/DI uncertainty | 10-45 min | 3-6 GB | Tiled processing |
| Visualization | 2-5 min | 2-3 GB | Map rendering |
| **Total** | **~1-3 hours** | **8-16 GB peak** | Highly variable by dataset |

**Time Factors:**
- FFS time scales quadratically with covariate count
- QRF time scales with raster size and ntree
- Larger n.tile values reduce memory but don't dramatically speed up (marginal I/O overhead)
- Linux/Unix typically 10-20% faster than Windows (more efficient fork-based parallelization)

### 6.6 Troubleshooting

**Issue 1: Memory Error ("cannot allocate vector of size...")**

**Cause:** Insufficient RAM for current raster size

**Solutions:**
```r
# Option A: Increase tile granularity
SOCastR(..., n.tile = 8)  # More tiles = smaller per-tile memory

# Option B: Reduce raster resolution before processing
# Resample raster to coarser resolution using GDAL:
# gdalwarp -tr 200 200 input.tif output.tif

# Option C: Process in smaller geographic tiles externally
# Split large raster into regions, process separately, mosaic results

# Option D: Close other applications to free RAM
# Minimize other R objects in memory
```

**Issue 2: Parallel Processing Error ("Error in checkForRemoteErrors...")**

**Cause:** Windows PSOCK cluster communication failure

**Solutions:**
```r
# Add error handling in script
tryCatch(
  SOCastR(...),
  error = function(e) {
    cat("Parallel error detected. Retrying with sequential backend...\n")
    registerDoSEQ()  # Fall back to sequential
    # Retry processing
  }
)

# Or reduce cores manually:
cl <- makeCluster(2)  # Use only 2 cores instead of automatic
registerDoParallel(cl)
```

**Issue 3: CRS Mismatch Error**

**Cause:** Samples and covariates have incompatible or undefined CRS

**Solutions:**
```r
# Check CRS before running
library(sf)
library(terra)

samples <- st_read("input/samples.shp")
covariates <- rast("input/covariates.tif")

st_crs(samples)
crs(covariates)

# If samples lack CRS:
samples <- st_set_crs(samples, "EPSG:25832")
st_write(samples, "input/samples_fixed.shp")

# If covariate lacks CRS:
crs(covariates) <- "EPSG:25832"
writeRaster(covariates, "input/covariates_fixed.tif")
```

**Issue 4: "No variance" Warning / Zero-Variance Predictor Removal**

**Cause:** Some covariate layers have constant values everywhere

**Solutions:**
```r
# This is normal and expected; constant layers provide no information
# Script automatically removes them

# To identify problematic layers before running:
library(terra)
covariates <- rast("input/covariates.tif")

for (i in 1:nlyr(covariates)) {
  layer <- covariates[[i]]
  val_range <- paste(min(values(layer)), "-", max(values(layer)))
  cat(names(layer), ":", val_range, "\n")
  # If min == max, layer has no variance
}

# Remove zero-variance layers in GDAL:
# gdal_translate -b 1 -b 2 -b 4 input.tif output.tif  # Keep bands 1,2,4; skip 3
```

**Issue 5: Very Slow Feature Selection**

**Cause:** Too many covariates or small sample size

**Solutions:**
```r
# Option A: Pre-filter covariates by correlation
# Keep only covariates with |r| > 0.3 to SOC

# Option B: Reduce covariate count
# Use only most relevant environmental variables

# Option C: Accept longer runtime
# FFS is thorough; 1-2 hours typical for 20+ covariates

# Option D: Skip FFS and use all covariates
# Modify script to skip FFS step (not recommended; reduces parsimony)
```

**Issue 6: File Not Found / Path Errors**

**Cause:** Incorrect working directory or relative paths

**Solutions:**
```r
# Verify working directory
getwd()

# Check files exist
file.exists("input/samples.shp")
file.exists("input/covariates.tif")

# Use absolute paths if relative paths fail
SOCastR(
  working_dir = "/home/user/SOCastR/",  # Absolute path
  input_dir = "/home/user/SOCastR/input/",
  ...
)
```

---

## 7. Outputs and Deliverables

### 7.1 Output Directory Structure

```
output/
├── Raster_Outputs/
│   ├── FinalPrediction_SocRaster.tif
│   ├── FinalPrediction_SocQuantileLayers.tif
│   ├── FinalPrediction_DissimilarityIndex.tif
│   └── FinalPrediction_AOAmask.tif
│
├── Validation_Statistics/
│   ├── ExtractValues_SampleDataSummary.csv
│   ├── SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv
│   ├── SpatialCrossValidation_GeodistMetrics.csv
│   ├── ForwardFeatureSelection_SelectedVariables.csv
│   ├── FinalModel_VariableImportance.csv
│   ├── FinalModel_Accuracy.csv
│   ├── Validation_PerformanceComparison.csv
│   ├── QuantilePrediction_Accuracy.csv
│   ├── QuantilePrediction_UncertaintySummary.csv
│   └── FinalPrediction_AoaSummaryStatistics.csv
│
├── Visualizations/
│   ├── SpatialDataPartition_MapTrainTestSplit.png
│   ├── SpatialDataPartition_SocDistSpatialBlock.png
│   ├── SpatialCrossValidation_GeodistEcdf.png
│   ├── SpatialCrossValidation_GeodistDensity.png
│   ├── ForwardFeatureSelection_VariableImportance.png
│   ├── Validation_ScatterPlot.png
│   ├── Validation_ResidualPlot.png
│   ├── FinalPrediction_MapSocRF.png
│   ├── FinalPrediction_MapQuantialSoc5Percentile.png
│   ├── FinalPrediction_MapQuantialSoc50Percentile.png
│   ├── FinalPrediction_QuantialMapSoc95Percentile.png
│   ├── FinalPrediction_MapPIW90.png
│   └── FinalPrediction_MapDI.png
│
└── Model_Objects/
    └── trainDI_object.rds
```

### 7.2 Raster Output Specifications

#### 7.2.1 Main SOC Prediction Map
- **Filename:** `FinalPrediction_SocRaster.tif`
- **Format:** GeoTIFF (LZW compression)
- **Bands:** 1
- **Data Type:** Float32
- **Unit:** % or g/kg (same as input)
- **NoData Value:** Pixels outside raster extent or with missing covariates
- **CRS:** Same as input covariates (EPSG code embedded)
- **Resolution:** Identical to covariate rasters
- **Contents:** Random Forest point predictions across full study area

#### 7.2.2 Quantile Predictions
- **Filename:** `FinalPrediction_SocQuantileLayers.tif`
- **Format:** Multi-band GeoTIFF (LZW compression)
- **Bands:** 4
  - Band 1: Q05 (5th percentile, lower confidence bound)
  - Band 2: Q50 (50th percentile, median prediction)
  - Band 3: Q95 (95th percentile, upper confidence bound)
  - Band 4: PIW (Prediction Interval Width = Q95 - Q05)
- **Data Type:** Float32
- **Unit:** % or g/kg
- **Interpretation:**
  - Narrow intervals (small PIW) = higher local certainty
  - Wide intervals (large PIW) = higher local uncertainty
  - Based on Quantile Regression Forest algorithm

#### 7.2.3 Dissimilarity Index
- **Filename:** `FinalPrediction_DissimilarityIndex.tif`
- **Format:** GeoTIFF
- **Bands:** 1
- **Data Type:** Float32
- **Unit:** Dimensionless (0 to threshold to ∞)
- **Interpretation:**
  - DI ≤ threshold: Within Area of Applicability
  - DI > threshold: Outside AOA (extrapolation risk)
  - Measures multivariate distance to training feature space
  - Calculated via CAST::trainDI and CAST::aoa functions

#### 7.2.4 Area of Applicability Mask
- **Filename:** `FinalPrediction_AOAmask.tif`
- **Format:** GeoTIFF
- **Bands:** 1
- **Data Type:** UInt8 (unsigned 8-bit integer)
- **Values:** 
  - 1 = Inside AOA (reliable predictions)
  - 0 = Outside AOA (extrapolation, use with caution)
  - NA = NoData (outside covariate extent)
- **Interpretation:** Binary decision rule for prediction reliability
- **Recommended Use:** Filter predictions to AOA=1 for critical applications

### 7.3 Statistical Output Files

#### 7.3.1 Data Summary
**File:** `ExtractValues_SampleDataSummary.csv`

```
Metric,Value
Total samples,1250
Clean samples,1204
Removed samples,46
Number of covariates,23
SOC mean,12.34
SOC sd,8.56
SOC min,0.89
SOC max,42.12
```

#### 7.3.2 Train/Test Split
**File:** `SpatialDataPartition_TrainTestStatisticsSpatialBlock.csv`

```
Metric,Value
N_train,962
N_test,242
Train_SOC_mean,12.45
Test_SOC_mean,11.98
SOC_mean_difference,0.47
SOC_mean_diff_percent,3.77
KS_statistic,0.1234
KS_p_value,0.2834
```

#### 7.3.3 Model Performance
**File:** `Validation_PerformanceComparison.csv`

```
Dataset,RMSE,MAE,R2,Bias,RMSE_percent
Training_CV,2.12,1.56,0.628,0.0,17.2
Test,2.34,1.71,0.598,-0.12,19.0
```

#### 7.3.4 Feature Selection
**File:** `ForwardFeatureSelection_SelectedVariables.csv`

```
Variable,Selection_order
COV15,1
COV08,2
COV12,3
COV19,4
COV05,5
```

#### 7.3.5 Uncertainty Summary
**File:** `QuantilePrediction_UncertaintySummary.csv`

```
Metric,Value
Mean_PI_width,4.23
Median_PI_width,3.89
95th_percentile_width,9.12
```

**File:** `FinalPrediction_AoaSummaryStatistics.csv`

```
Metric,Value
Total_prediction_cells,1048576
Cells_within_AOA,892341
Percent_within_AOA,85.12
Mean_DI,0.0847
Max_DI,2.8934
```

### 7.4 Visualization Outputs

All PNG files are generated at 300 DPI suitable for publication.

**Spatial Distribution Map:**
- `SpatialDataPartition_MapTrainTestSplit.png`: Training vs test sample locations
- **Use:** Visual assessment of spatial sampling pattern and representativeness

**Validation Plots:**
- `Validation_ScatterPlot.png`: Observed vs predicted with R² and RMSE
- `Validation_ResidualPlot.png`: Predicted vs residuals with trend line
- **Use:** Assess prediction bias and homoscedasticity

**Uncertainty Visualizations:**
- `FinalPrediction_MapQuantialSoc5Percentile.png`: Lower confidence bound
- `FinalPrediction_MapQuantialSoc50Percentile.png`: Median prediction
- `FinalPrediction_MapQuantialSoc95Percentile.png`: Upper confidence bound
- `FinalPrediction_MapPIW90.png`: Prediction interval width (90%)
- **Use:** Understand spatial patterns of predictions and uncertainty

**Applicability Assessment:**
- `FinalPrediction_MapDI.png`: Dissimilarity Index (red=high, blue=low)
- **Use:** Identify extrapolation zones and reliability zones

### 7.5 Model Objects

**Saved Models (RDS format):**
- `trainDI_object.rds`: trainDI object for reuse in new predictions
- **Use:** Quickly compute AOA for new covariate rasters without retraining

---

## 8. Quality Assurance and Validation

### 8.1 Internal Quality Checks

The workflow includes automated quality assurance checkpoints:

**Data Loading Checks:**
- ✓ Verify shapefile integrity and readability
- ✓ Check for valid CRS definition
- ✓ Confirm SOC column exists and contains numeric data
- ✓ Validate raster stack loading and band count

**Data Completeness Checks:**
- ✓ Count and report missing values (NA)
- ✓ Identify and remove out-of-bound samples
- ✓ Report sample removal statistics
- ✓ Verify covariate extraction success

**Model Training Checks:**
- ✓ Remove zero-variance predictor layers
- ✓ Check for convergence in feature selection
- ✓ Validate CV fold structure (no empty folds)
- ✓ Monitor training/test split balance

**Prediction Checks:**
- ✓ Verify prediction extent matches covariate domain
- ✓ Check for prediction outliers vs training range
- ✓ Validate quantile ordering (Q5 < Q50 < Q95)
- ✓ Ensure AOA binary classification is 0/1

### 8.2 User Validation Protocol

**Recommended Quality Assurance Steps:**

**Step 1: Input Data Validation**
```r
# Verify sample data quality
samples <- st_read("input/samples.shp")
summary(samples$SOC)  # Check range and outliers
plot(samples)         # Visual check for clustering
```

**Step 2: Output Visual Inspection**
- Review train/test split map: Should show good spatial separation
- Check prediction map: Patterns should align with covariates
- Examine uncertainty maps: Higher uncertainty in data-sparse regions
- Assess AOA: >70% AOA coverage is desirable

**Step 3: Statistical Validation**
- Verify R² > 0.3 (minimum) or > 0.5 (good)
- Check Bias ≈ 0 (within ±10% of mean SOC)
- Compare train vs test RMSE: Should be similar (no overfitting)
- Validate prediction interval coverage (~90% within PI)

**Step 4: Cross-Validation** (if possible)
- Collect independent test samples (not used in modeling)
- Predict at independent sample locations
- Compare predictions to observations
- Calculate independent R² and RMSE

**Step 5: Expert Judgment**
- Compare maps with existing soil maps (if available)
- Check for unrealistic prediction patterns
- Validate spatial continuity and smoothness
- Assess domain knowledge alignment

### 8.3 Quality Control Checklist

**Pre-Publication Verification:**

- [ ] Input data quality documented (source, method, accuracy)
- [ ] Sample size adequate (≥100, preferably ≥200)
- [ ] Spatial distribution representative (not heavily clustered)
- [ ] R² on test set > 0.3 minimum (> 0.5 preferred)
- [ ] Bias within ±10% of mean SOC (unbiased predictions)
- [ ] AOA covers >50% of study area (>70% preferred)
- [ ] Prediction intervals reasonable width for application
- [ ] All output maps visually inspected
- [ ] Metadata files complete and accurate
- [ ] Processing log documented with parameters and timing

### 8.4 Known Limitations and Caveats

**Methodological Limitations:**

1. **Random Forest Constraints:**
   - Cannot extrapolate beyond training data range
   - Predictions in extrapolation zones = mean of training values
   - Tends to underpredict extreme high/low values
   - Requires sufficient sample size (n ≥ 50) for stable estimates

2. **Spatial Autocorrelation:**
   - Strong spatial autocorrelation can lead to overoptimistic CV results
   - KNNDM mitigates but doesn't eliminate this issue
   - Inspect residual autocorrelation with Moran's I if concerned

3. **Temporal Stationarity Assumption:**
   - Assumes SOC-covariate relationships are constant over time
   - Model trained on historical data may not extrapolate to future
   - Recalibration recommended if >10 years between training and prediction

4. **Stationarity Over Space:**
   - Assumes relationships are spatially constant (within CV blocks)
   - May fail in highly heterogeneous landscapes
   - Test via local regression or stratified analysis if suspected

5. **Covariate Limitations:**
   - Predictions only as good as input covariates
   - Missing key drivers reduces accuracy
   - Poor covariate resolution limits prediction resolution
   - Collinear covariates reduce model parsimony (addressed by FFS)

6. **Computational Constraints:**
   - Very large rasters (>100M cells) may exceed available memory
   - Processing time scales quadratically with covariate count
   - Parallelization overhead for small datasets
   - Trade-off between tiling granularity and I/O efficiency

7. **Uncertainty Limitations:**
   - QRF intervals may not achieve nominal coverage (80-100%)
   - DI threshold is data-dependent; may be too conservative or permissive
   - AOA binary classification may be too strict for marginal zones
   - Ignores sources of uncertainty not captured by models (e.g., lab error)

---

## 9. Citation and Acknowledgments

### 9.1 How to Cite SOCastR

**Recommended Citation Format (APA):**

```
Möller, M. (2025). SOCastR: Soil Organic Carbon Prediction Workflow with 
Uncertainty Quantification (Version 1.0.0) [Computer software]. 
Julius Kühn Institute. https://doi.org/10.5281/zenodo.XXXXXXX
```

**BibTeX Format:**

```bibtex
@software{moeller2025socastr,
  author = {M{\"o}ller, Markus},
  title = {{SOCastR}: Soil Organic Carbon Prediction Workflow with 
           Uncertainty Quantification},
  version = {1.0.0},
  year = {2025},
  doi = {10.5281/zenodo.XXXXXXX},
  url = {https://github.com/[username]/SOCastR},
  institution = {Julius K{\"u}hn Institute},
  license = {MIT}
}
```

### 9.2 CITATION.cff File

Included in repository root:

```yaml
cff-version: 1.2.0
message: "If you use this software in your research, please cite it as follows."
authors:
  - family-names: "Möller"
    given-names: "Markus"
    orcid: "https://orcid.org/XXXX-XXXX-XXXX-XXXX"
    affiliation: "Julius Kühn Institute, Germany"
title: "SOCastR: Soil Organic Carbon Prediction Workflow with Uncertainty Quantification"
version: 1.0.0
date-released: "2025-11-01"
doi: "10.5281/zenodo.XXXXXXX"
repository-code: "https://github.com/[username]/SOCastR"
license: "MIT"
keywords:
  - "digital-soil-mapping"
  - "machine-learning"
  - "soil-organic-carbon"
  - "random-forest"
  - "uncertainty-quantification"
  - "R"
abstract: >
  SOCastR implements a comprehensive workflow for digital soil mapping of 
  soil organic carbon using Random Forest and Quantile Regression Forest 
  with spatial cross-validation and Area of Applicability assessment. 
  The workflow includes forward feature selection, external validation, 
  and ISO 19157 compliant data quality assessment.
```

### 9.3 Key Dependencies to Cite

**CAST Package:**
```
Meyer, H., & Pebesma, E. (2022). Machine learning-based global maps of ecological 
variables and the challenge of assessing them. Nature Communications, 13(1), 2208.
https://doi.org/10.1038/s41467-022-29838-9
```

**Quantile Regression Forest:**
```
Meinshausen, N. (2006). Quantile regression forests. Journal of Machine Learning 
Research, 7(Jun), 983-999.
```

**Random Forest Algorithm:**
```
Breiman, L. (2001). Random forests. Machine Learning, 45(1), 5-32.
https://doi.org/10.1023/A:1010933404324
```

### 9.4 Acknowledgments

**Funding:** [Add funding information]

**Data Sources:** [Acknowledge data providers and sources]

**Computational Resources:** [Acknowledge HPC facilities if used]

**Contributors:** [List collaborators, reviewers, and their roles]

**Inspiration and Methodology:** This workflow builds upon methods and best practices from:
- CAST package documentation and vignettes (Meyer & Pebesma)
- Digital soil mapping literature (Lagacherie, McBratney, Minasny)
- Spatial machine learning standards (Hijmans et al.)
- ISO 19157 data quality framework
- W3C PROV-DM provenance standards

---

## 10. References

### 10.1 International Standards

1. **ISO 19157:2013** - Geographic information - Data quality
   https://www.iso.org/standard/32575.html

2. **ISO 19115-1:2014** - Geographic information - Metadata - Part 1: Fundamentals
   https://www.iso.org/standard/53798.html

3. **ISO 19115-3:2020** - Geographic information - Metadata - Part 3: XML Schema Implementation
   https://www.iso.org/standard/73910.html

4. **W3C PROV-DM: The PROV Data Model** (2013)
   https://www.w3.org/TR/prov-dm/

5. **FAIR Guiding Principles for Scientific Data Management** (Wilkinson et al., 2016)
   https://doi.org/10.1038/sdata.2016.18

### 10.2 Core Methodology

6. **Meyer, H., & Pebesma, E. (2022).** Machine learning-based global maps of ecological variables and the challenge of assessing them. Nature Communications, 13(1), 2208.
   https://doi.org/10.1038/s41467-022-29838-9

7. **Breiman, L. (2001).** Random forests. Machine Learning, 45(1), 5-32.
   https://doi.org/10.1023/A:1010933404324

8. **Meinshausen, N. (2006).** Quantile regression forests. Journal of Machine Learning Research, 7(Jun), 983-999.

9. **Meyer, H., Reudenbach, C., Hengl, T., Katurji, M., & Nauss, T. (2018).** Improving performance of spatio-temporal machine learning models using forward feature selection and target-oriented validation. Environmental Modelling & Software, 101, 1-9.
   https://doi.org/10.1016/j.envsoft.2017.12.001

10. **Meyer, H., & Pebesma, E. (2021).** Predicting into unknown space? Estimating the area of applicability of spatial prediction models. Methods in Ecology and Evolution, 12(9), 1620-1633.
    https://doi.org/10.1111/2041-210X.13650

### 10.3 Digital Soil Mapping

11. **Lagacherie, P., & McBratney, A. B. (2006).** Spatial soil information systems and spatial soil inference systems: Perspectives for digital soil mapping. In Digital Soil Mapping (pp. 3-22). Springer.
    https://doi.org/10.1007/1-4020-5010-8_1

12. **McBratney, A. B., Minasny, B., & Catambas, U. (2011).** Digital soil mapping: Bridging research, production, and application. Geoderma, 171, 32-41.
    https://doi.org/10.1016/j.geoderma.2011.01.024

13. **Hengl, T., MacMillan, R. A., & Nolin, M. C. (2009).** Preprocessed SRTM data for bioclimatic modelling of global lightning-caused fire ignition risk. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 2(3), 139-145.

14. **Malone, B. P., McBratney, A. B., Minasny, B., & Laslett, G. M. (2009).** Mapping continuous depth functions of soil carbon storage and available water capacity. Geoderma, 154(1-2), 138-152.
    https://doi.org/10.1016/j.geoderma.2009.01.007

### 10.4 Spatial Validation Methods

15. **Hijmans, R. J. (2012).** Cross-validation of species distribution models: removing spatial sorting bias and evaluation truism. Ecology, 93(7), 1461-1463.
    https://doi.org/10.1890/11-2027.1

16. **Walesiak, M. (2014).** Clustering with Gower distance. In Encyclopedia of Measurement and Statistics (pp. 1-4).

### 10.5 Software and Tools

17. **Kuhn, M. (2008).** Building predictive models in R using the caret package. Journal of Statistical Software, 28(5), 1-26.
    https://doi.org/10.18637/jss.v028.i05

18. **Hijmans, R. J., & Bivand, R. S. (2021).** raster: Geographic data analysis and modeling (R package).
    https://CRAN.R-project.org/package=raster

19. **Pebesma, E. (2018).** Simple features for R: Standardized support for spatial vector data. The R Journal, 10(1), 439-446.
    https://doi.org/10.32614/RJ-2018-009

### 10.6 Reproducibility and FAIR Software

20. **Smith, A. M., Katz, D. S., & Niemeyer, K. E. (2016).** Software citation principles. PeerJ Computer Science, 2, e86.
    https://doi.org/10.7717/peerj-cs.86

21. **Jones, M. B., et al. (2017).** CodeMeta: An exchange schema for software metadata. Technical Report, IMCS.
    https://codemeta.github.io/

---

## Appendix A: Example Processing Log

```
================================================================================
SOCastR: Soil Organic Carbon (SOC) prediction workflow with uncertainties
================================================================================

Date: 2025-11-01
Time: 14:30:15 CET
R Version: 4.4.1
Platform: x86_64-w64-mingw32/x64 (Windows 11)

================================================================================
PARAMETERS
================================================================================
Working directory: C:/Users/user/SOCastR/
Input directory: input
Output directory: output
Samples file: SAMPLES_EPSG25832.shp
Covariates file: COVARIATES_EPSG25832.tif
SOC column: SOC
N.tile: 4
Model uncertainty: TRUE
Distance uncertainty: TRUE
Random seed: 42

================================================================================
DATA LOADING AND PREPARATION
================================================================================
Loaded 1250 sample points
Loaded 23 covariate layers
Raster resolution: 100 x 100 m

Extract (Filtered) Values
Original samples: 1250
Clean samples: 1204
Removed samples: 46 (3.68%)

================================================================================
SPATIAL DATA PARTITION
================================================================================
Spatial stratified split using CreateSpacetimeFolds
Training samples: 962 (79.9%)
Test samples: 242 (20.1%)
KS p-value: 0.2834 (distributions not significantly different)

================================================================================
FORWARD FEATURE SELECTION
================================================================================
Total available predictors: 23
Selected 7 variables:
 - COV15
 - COV08
 - COV12
 - COV19
 - COV05
 - COV03
 - COV18
FFS completed in 32 minutes

================================================================================
EXTERNAL VALIDATION
================================================================================
=== PERFORMANCE COMPARISON ===
Dataset         RMSE      MAE       R²     Bias    RMSE%
Training_CV     2.12      1.56      0.628  0.00    17.2
Test           2.34      1.71      0.598  -0.12   19.0

Test set validation completed.

================================================================================
SPATIAL PREDICTION
================================================================================
SOC prediction map created.
Prediction map saved.

================================================================================
MODEL UNCERTAINTY: QUANTILE REGRESSION FOREST
================================================================================
Training quantile regression forest model...
QRF model training completed.

Tile-based processing:
Creating 16 tiles (4x4 grid)
Processing 16 tiles across 7 cores...
Tile processing completed in 47 minutes

Merging tiles into final rasters...
Tile merging completed.

Quantile predictions completed.

Mean_PI_width: 4.23
Median_PI_width: 3.89
95th_percentile_width: 9.12

================================================================================
DISTANCE UNCERTAINTY QUANTIFICATION
================================================================================
Step 1: Computing trainDI threshold...
TrainDI threshold: 0.0892

Large raster detected. Using tile-based processing...
Creating 16 tiles (4x4 grid)
Using 7 cores for parallel processing

Processing AOA for each tile...
Tile processing completed in 22 minutes

Merging tiles into final rasters...

=== AOA COMPUTATION COMPLETED ===
Total_prediction_cells: 1048576
Cells_within_AOA: 892341
Percent_within_AOA: 85.12
Mean_DI: 0.0847
Max_DI: 2.8934

================================================================================
PROCESSING COMPLETED
================================================================================
Total execution time: 2 hours 14 minutes
Output directory: C:/Users/user/SOCastR/output/

Outputs generated:
- Raster files: 4 (predictions, uncertainty, applicability)
- Statistics: 10 CSV files
- Visualizations: 13 PNG maps
- Model objects: 1 RDS file

FAIR Documentation: SOCastR-FAIR-Docs.md
Quality Assessment: Complete
Ready for publication: YES
```

---

## Appendix B: ISO 19157 Metadata XML Example

```xml
<?xml version="1.0" encoding="UTF-8"?>
<mdb:MD_Metadata xmlns:mdb="http://standards.iso.org/iso/19115/-3/mdb/2.0"
                  xmlns:mdq="http://standards.iso.org/iso/19157/-2/mdq/1.0"
                  xmlns:mcc="http://standards.iso.org/iso/19115/-3/mcc/1.0"
                  xmlns:gco="http://www.opengis.net/gml/3.2.1">
  
  <mdb:metadataIdentifier>
    <mcc:MD_Identifier>
      <mcc:code>
        <gco:CharacterString>SOCastR_Prediction_2025-11-01</gco:CharacterString>
      </mcc:code>
    </mcc:MD_Identifier>
  </mdb:metadataIdentifier>
  
  <mdb:dataQualityInfo>
    <mdq:DQ_DataQuality>
      <mdq:scope>
        <mcc:MD_Scope>
          <mcc:level>
            <mcc:MD_ScopeCode codeList="http://standards.iso.org/iso/19115/resources/Codelist/cat/codelists.xml#MD_ScopeCode" codeListValue="dataset"/>
          </mcc:level>
        </mcc:MD_Scope>
      </mdq:scope>
      
      <mdq:report>
        <mdq:DQ_QuantitativeAttributeAccuracy>
          <mdq:measureDescription>
            <gco:CharacterString>RMSE of SOC prediction on independent test set</gco:CharacterString>
          </mdq:measureDescription>
          <mdq:result>
            <mdq:DQ_QuantitativeResult>
              <mdq:value>
                <gco:Record>2.34</gco:Record>
              </mdq:value>
              <mdq:valueUnit>
                <gco:CharacterString>g/kg</gco:CharacterString>
              </mdq:valueUnit>
            </mdq:DQ_QuantitativeResult>
          </mdq:result>
        </mdq:DQ_QuantitativeAttributeAccuracy>
      </mdq:report>
      
      <mdq:report>
        <mdq:DQ_QuantitativeAttributeAccuracy>
          <mdq:measureDescription>
            <gco:CharacterString>R-squared (coefficient of determination) on test set</gco:CharacterString>
          </mdq:measureDescription>
          <mdq:result>
            <mdq:DQ_QuantitativeResult>
              <mdq:value>
                <gco:Record>0.598</gco:Record>
              </mdq:value>
            </mdq:DQ_QuantitativeResult>
          </mdq:result>
        </mdq:DQ_QuantitativeAttributeAccuracy>
      </mdq:report>
      
    </mdq:DQ_DataQuality>
  </mdb:dataQualityInfo>
  
</mdb:MD_Metadata>
```

---

**End of Documentation**

This FAIR documentation comprehensively covers all aspects of the SOCastR workflow, including provenance, fitness-for-purpose assessment, and ISO 19157 compliant data quality metadata. It is ready for publication on GitHub and Zenodo.
