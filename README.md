# Biological Age Project for Python README

**Author:** Elisa Warner  
**Date of Last Update:** March 29, 2020  
**Contact: elisawa@umich.edu**  
**Version: 2.0**  
**Note: If you want any customizations to the code, I am willing to do this for free with credit to the author**

## Description
This project is based off the paper "Modeling the Rate of Senescence: Can Estimated Biological Age Predict Mortality More Accurately Than Chronological Age?" by Morgan Levine (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3660119/), which is based on the Klemera-Doubal Method for calculating biological age (https://www.ncbi.nlm.nih.gov/pubmed/16318865). Featured below is an implementation in Python. The original paper was written to stratify males and females, but it is recognized here that this may not be the case for all implementations. Therefore, you can create "One Group" version is available, which includes no stratification (See _How to Use_).

## New In This Version:
Version 2.0 contains the following updates from the previous version:
* Removed SAS implementation
* Removed Python 2.7 implementation
* Included automatic group detection for up to 2 groups
* Allow inputs to be DataFrames instead of a list
* Included KDM_model as built-in function for the function library, rather than a function you have to write
* Included a result_stat function which summarizes results against age
* Removed ability of graph to plot SLR of age against age, seqn, sex, samp_wt

## What to Expect
1. **Input:** 
* A DataFrame of subjects, where the features for each subject is thought to predict age  
* The name you would like to assign the cache. Or the name of an existing cache that contains the parameters for your data
* The index of the "age" column  
* The index of the "sample weight" column  
* The index of the "sex" or "group" column (note that the code assumes 1 or 2 groups)  
* The index of the primary key column (a subject identifier)  
* The name you want to assign to the output file
2. **Output:** 
* A saved dataset file of the original data and two extra columns: 1) the calculated biological age of the subject, 2) the corrected biological age of the subject, utilizing the true chronological age of the subject.  
* A cache file that contains the parameters of all the regressions needed  
* Results [list]: residual stats of corrected biological age (BAC) compared against the true age of the individual  

## Background
In the original study, Levine applies the NHANES III dataset to the Klemera-Doubal Method. She first calculates a biological age of an individual based on the following equation:
![](images/BA.png)
where $q$ represents the intercept, $k$ represents the slope, and $s$ represents the MSE of every feature column $j$ regressed against age.  

She then calculates the Corrected Biological Age of the same individual based on the following equation:
![](images/BAC.png)
where $CA$ represents the chronological age (true age) of the individual, and $s^2$ represents an age- and r-value-corrected MSE.

## Files Included:
1. bio_age.ipynb : Python Notebook version
2. bio_age.py : python script version
4. BioAgeSample_NHANESIII_without_fvc.txt : NHANES III sample dataset give to me by the author.
5. Example.ipynb : A python notebook that gives an example of how to use the code.

## Requirements
1. To run bio_age.py : Python 3.6+ with Anaconda
2. To run bio_age.ipynb : Python 3.6+ with Anaconda and Jupyter Lab

**You will need to install the following packages:**
* numpy  
* math  
* scipy  
* scipy.stats  
* sklearn  
* pandas  
* json  
* matplotlib.pyplot  

## Important Notes: 
1. Note that the cache is intended to save the calculated parameters. Therefore, if you change your dataset and need to update parameters, you must DELETE the cache file(s).  
2. For best use of biological age in general, it is important to run correlations on your dataset features to include only those features which are significantly correlated with age.

## How to Use:
### Dataset Setup (INPUT)
The original paper utilizes sample weights and stratifies by sex of the individual. Therefore, the dataset INPUT requires four special data columns:  
1. `seqn` : the unique primary key of each observation in the dataset  
2. `group` : this can be the sex of the subject (binary 0 or 1) or all '1's if you don't want stratification  
3. `samp_wt` : these are the sample weights for each individual (range: [0 1]). If you don't want to include sample weights, set this column with all '1's.  
4. `age` : This column is required for this study.  

Note that columns do not have to be named this, but you must identify in the input code to which column index of the DataFrame each of these columns belongs. To use this code, the program requires a delimited text file (either .csv, .tsv, or .txt). The current version requires datasets to be in DataFrame form.

### Example Code 1
In the most simple example, our dataset `df` is a DataFrame set up according to _Dataset Setup_ above. And we have one dataset (not split into training and testing). The variables `age_index`, `genderindex`, `seqn`, `samp_wt_index` are the python indices of the age, group, seqn, and sample_wt columns, respectively. 

**Cache**  
`cachename` is the name of a cache file. If no cache file exists yet, assign any name. The cache file will be created automatically. It saves all the calculated parameters so that training and testing can be possible (the biological ages of the test set are simply based off of the saved parameters from the training set's cache file). The `cachename` argument allows users to save multiple cache files for more efficient execution.

**Execution Code:**

    import bio_age as BA  # import script

    model = BA.Methods(df, cachename, age = age_index, genderindex = genderindex, primarykey = seqn, samp_wt = samp_wt_index, 'Training_Results.csv')`
    results = model.KDM() # run model`

### Example Code 2
In this example, we run the built-in function for the KDM model. The assumption is that we have a training set (`trainset`) and testing set (`testset`), and an assigned Cache File (doesn't have to exist yet) (`cachename`). The training set and testing set are both pandas DataFrames.

The output gives some statistics of how the model performs against the testing set.

    # find indices for age, group, primary index, sample weights
    header_names = list(trainset)
    age_index = header_names.index('age')
    genderindex = header_names.index('GROUP') # this is the name of the group
    seqn = header_names.index('seqn')
    samp_wt_index = header_names.index('samp_wt')

    results = BA.KDM_model(trainset, trainset, 'cache_ex2.json','MouseAge_results', age_index, genderindex, seqn, samp_wt_index,'Output_Ex1')

** Note: When you're using the KDM_model() function, the output filename should just be a stem and NOT contain ('.csv') at the end. The reason for this is that '\_train.csv' and '\_test.csv' are added to the end of the output filename.

## Output Explained
An example output, saved to a csv file, will appear like so:
![results](images/Output_ex.png)
The complete dataset is still visible, but two columns appear at the end of the file: `BA` and `BAC`. `BA` is defined above in _Background_ as the initial, uncorrected Biological Age. `BAC` is also defined above as the corrected Biological Age based both on `BA` and on the true chronological age of the individual.

## Hyperparameters
There are two hyperparameters in this code:  
1. `GRAPHON` : With `GRAPHON = True`, the simple linear regression plots for each feature regressed against age will show upon calculation of the linear regression (calculations only occur when the cache file does not exist). This option is turned on by default (`GRAPHON = True`), but they won't display if parameters are loaded from the cache.
2. `AGEON` : This is a linear transformation of the SBA variable (_see paper_) that was conducted in the original paper  
            # Extra step: Linearly transform so that SBA maintains same mean but now linearly increases with age, so difference is 5 between CAmax and CAmin
In my experience, this transformation does little to change the overall results. It is set by default to true (`AGEON = True`)  

These hyperparameters are hard to find in the notebook code. However, to change them, find the following lines, and change the values as needed. Note that if `GRAPHON = True`, `import ggplot as gg` also needs to be uncommented:

        self.AGEON = True # age corrector toggle
        self.GRAPHON = False # do you want to see linear regression graphs for each variable on age?

## Extra Features
There are some additional features that can be viewed with this code (calculate and display all means, correlations with age, and view results in graphical form). View examples of how to use these by viewing the Example.ipynb notebook.