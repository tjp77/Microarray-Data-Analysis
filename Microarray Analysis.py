# Part III, KNN
# https://scikit-learn.org/stable/modules/generated/sklearn.svm.libsvm.fit.html#sklearn.svm.libsvm.fit
def ClassifyGenes(genesTraining, labelsTraining, genesTesting, labelsTesting):
    
    knnSize = 0.6 * len(genes);
    
    knn = KNeighborsClassifier(algorithm = 'auto', leaf_size = 30, # what was leaf size? 
                               metric = 'minkowski', metric_params = None, n_jobs = 1,
                               n_neighbors = 3, p = 2, weights = 'uniform');
    
    knn.fit(genesTraining, labelsTraining);
    
    print("Classifier predictions for testing data:");
    print(knn.predict(genesTesting));
    
    print("Target labels for testing data:"); # --- Save output to file for easier read 
    print(labelsTesting);
    
    return 0;


def WriteClassifyResultsToFile(testingPrediction, labelsTesting):
    
    fileName = "classificationResults.txt";
        
    try:
        
        file = open(fileName, "w");
            
        for i in range(0, len(testingPrediction)):
                
            file.write(testingPrediction + "\t");
        
        file.write("\n");
        
        for i in range(0, len(testingPrediction)):
                
            file.write(testingPrediction + "\t");
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nClassification result save successful.\n");
            
    file.close();
    return 0;


def main():
    
    inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
    genesTraining = inputTraining[0];
    typesTraining = inputTraining[1];
    #print(types);
    
    Preprocess(genesTraining); 
    print("\nProcessed.\n");
    SaveAffymetrics1(genesTraining, typesTraining);
    
    inputTesting = readInData("Loading Leuk_ALL_AML.test .....", "Leuk_ALL_AML.test.txt", 20);
    genesTesting = inputTesting[0];
    typesTesting = inputTesting[1];
    
    SaveAffymetrics3(genesTesting, typesTesting);
    
    
    ClassifyGenes(genesTraining, typesTraining, genesTesting, typesTesting)
    
    return 0;
    


main();


# ======= To Do ======= 

# PreProcess function:
# Eliminate the genes with less than two fold change across the experiments (max/min <2);
# - REVIEW SAVED AFFY FILE, SEE IF SEEMS LIKE RIGHT THINGS REMOVED. ('seems' due to file too large to check all)


# Sort genes by p-values/T test stuff Part II 4. a

# output classifygenes prints to file for easier viewing and review.  - Needs testing

# have testing genes, need to, after picking the top 50 from Training, take just though matching genes (second excel file colums matching values) and use those only in the KNN part. [!!!] also only send those 50 to the SaveAffymetrics3 file.

# Complete KNN stuff.


# Once everything else is done, try dif amounts of n_neighbors for the KNN classifier and compare results. 



# ======= Finished From To Do =======

# - Output read in file and print to test if coming in ok.
# - Test/review SaveAffymetrics1 produced data compared to project example out put shown. 
# - - Opened in excel, format matches expected fine. 
# - - Types being filled ok. 




# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#
