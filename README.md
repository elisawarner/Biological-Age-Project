# Biological Age Project README

**Author:** Elisa Warner  
**Date of Last Update:** Dec 23, 2019  

**Description:** This project is based off the paper "Modeling the Rate of Senescence: Can Estimated Biological Age Predict Mortality More Accurately Than Chronological Age?" by Morgan Levine (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3660119/), which is based on the Klemera-Doubal Method for calculating biological age (https://www.ncbi.nlm.nih.gov/pubmed/16318865). It features my original SAS implementation of the project, then a cleaner python script. The original paper was written to stratify males and females, but it is recognized here that this may not be the case for all implementations. Therefore, a "One Group" version is available, which includes no stratification.

## Background
In the original study, Levine applies the NHANES III dataset to the Klemera-Doubal Method. She first calculates a biological age of an individual based on the following equation:
![](BA.png)
where $q$ represents the intercept, $k$ represents the slope, and $s$ represents the MSE of every feature column $j$ regressed against age.  

She then calculates the Corrected Biological Age of the same individual based on the following equation:
![](BAC.png)
where $CA$ represents the chronological age (true age) of the individual, and $s^2$ represents an age- and r-value-corrected MSE.

## Files Included:
1. BA_NB_Final_V2.ipynb
2. StudyReplicationReversedF.sas
3. BioAgeSample_NHANESIII_without_fvc.txt

## How to Use:
### Dataset Setup (INPUT)
The original paper utilizes sample weights and stratifies by sex of the individual. Therefore, the dataset INPUT requires four special data columns:
1. 'seqn' : the unique primary key of each observation in the dataset  
2. 'group' : this can be the sex of the subject (binary 0 or 1) or all '1's if you don't want stratification  
3. 'samp_wt' : these are the sample weights for each individual (range: [0 1]). If you don't want to include sample weights, set this column with all '1's.  
4. 'age' : This column is required for this study.  

To use this code, the program requires a delimited text file (either .csv, .tsv, or .txt). The current version requires datasets to be in list form rather than DataFrame. "List form" means that the data must be arranged as a list of lists, where the first row of the list of lists is the list of column names. If you start with a dataframe, this can be easily created by using Pandas command `to_list()`:  

**Convert DataFrame to list form**  

    # change format of train set to old format that BA method requires
    header_names = list(trainset)
    trainset_list = trainset.values.tolist()
    trainset_list.insert(0, header_names)

### Example Code 1
In the most simple example, our dataset is set up according to **Dataset Setup** above. And we have one dataset (not split into training and testing). The variables `age_index`, `genderindex`, `primaryindex`, `samp_wt_index` are the python indices of the age, group, seqn, and sample_wt columns, respectively. 

**Cache**  
`cachename` is the name of a cache file. If no cache file exists yet, assign any name. The cache file will be created automatically. It saves all the calculated parameters so that training and testing can be possible (the biological ages of the test set are simply based off of the saved parameters from the training set's cache file). The `cachename` argument allows users to save multiple cache files for more efficient execution.

**Execution Code:**

    import BA_NB_Final_V2 as BA`
    model = BA.Methods(dataset_list, cachename, age = age_index, genderindex = genderindex, primarykey = primaryindex, samp_wt = samp_wt_index)`
    results = model.KDM() # run model`

### Example Code 2
In this example, we create a function for the KDM model. The assumption is that we have a training set (`trainset`) and testing set (`testset`), and an assigned Cache File (doesn't have to exist yet) (`cachename`). The training set and testing set are both pandas DataFrames.

The output is simply a float array of all the corrected Biological Ages for each subject.


    def KDM_model(trainset, testset, cachename):
        # set up dataset columns
        trainset['seqn'] = range(1,trainset.shape[0]+1)
        trainset['group'] = 1
        trainset['samp_wt'] = 1

        testset['seqn'] = range(1,testset.shape[0]+1)
        testset['group'] = 1
        testset['samp_wt'] = 1
    
        testset = testset.drop(['GROUP', 'Days_alive', 'age_months'], axis = 1)
        trainset = trainset.drop(['GROUP', 'Days_alive', 'age_months'], axis = 1)
        #print(testset.shape, trainset.shape)
    
        # convert DataFrame to List Format
        header_names = list(trainset)
        trainset_list = trainset.values.tolist()
        trainset_list.insert(0, header_names)
    
        testset_list = testset.values.tolist()
        testset_list.insert(0, header_names)
    
        # find indices for age, group, primary index, sample weights
        age_index = header_names.index('age_days')
        genderindex = header_names.index('group')
        primaryindex = header_names.index('seqn')
        samp_wt_index = header_names.index('samp_wt')
    
        # Build model
        model = BA.Methods(trainset_list, cachename, age = age_index, genderindex = genderindex, primarykey = primaryindex, samp_wt = samp_wt_index)
        model.KDM() # train model
    
        # change format of test set from DataFrame to List Format
        header_names = list(testset)
        testset_list = testset.values.tolist()
        testset_list.insert(0, header_names)
    
        # test using trained model
        test_model = BA.Methods(testset_list, cachename, age = age_index, genderindex = genderindex, primarykey = primaryindex, samp_wt = samp_wt_index)
        results = test_model.KDM()
        results = pd.DataFrame(results[1:], columns = results[0])
        # test_results = test_results.append(results) # this was designed for 5-fold cross validation, so it would append the results of each fold into one dataframe

        # calculate stats
        stats = return_stats(results['BAC'].astype('float64'), results['age_days'].astype('float64')) # output Biological Age, Corrected
    
    return (stats)

### Output Explained

### Hyperparameters
