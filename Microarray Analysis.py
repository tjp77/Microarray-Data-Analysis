#!/user/bin/python
import sys;
import numpy as np
from sklearn import datasets
from sklearn.neighbors import KNeighborsClassifier

class Gene:
    
    def __init__(self, acc, all, aml, rowID):
        
        self.accessionId = acc;
        self.allList = all;
        self.amlList = aml;
        self.ID = rowID;


def readInData(prompt, fileName, threshold):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        try:
            
            print(prompt); 
            
            with open(fileName, "r") as file:
                
                firstLine = file.readline().split("\t"); 
                
                #for i in range(0, len(firstLine)):
                 #   print (str(i) + firstLine[i] + "|");
                
                # Not sure will need these, but have to pass through them either way.  
                line = file.readline(); 
                line = file.readline(); 
                
                # For first pass of while loop. 
                line = file.readline();
                # Keep track of which file line gene is on, so can easily pull testing genes corresponding to selected top 50 training genes. 
                id = 0;
                
                # Go through remaining lines and get all things interested in.
                while line:
                    
                    line = line.split("\t")
                    
                     # For all where first col, description contains "(endogenous control)", skip since supossed to eliminate them.
                    if ("endogenous control" in line[0]):
                        
                        line = file.readline();
                        continue;
                        
                    # second column is Accession ID.
                    accession = line[1];
                    all = [];
                    aml = [];
                    
                    j = 2;
                    # Get each expresson value.
                    while (j < len(line) - 1):
                        
                        if ( (int)(line[j]) < threshold ):
                            line[j] = threshold; 
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            all.append([(int)(line[j]), line[j+1]]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            aml.append([(int)(line[j]), line[j+1]]);
                        
                        j += 2;
                    
                    genes.append(Gene(accession, all, aml, id)); #print(genes[-1].allList, "\n\n");
                    id += 1;
                    line = file.readline(); 
        except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            sys.exit();
            
        else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return [genes, firstLine];

# Save Training Data, Part I
def SaveAffymetrics1(genes, types):
    
    fileName = "affymetrics1.txt";
        
    try:
            
        file = open(fileName, "w");
            
        expCount = len(genes[0].allList) + len(genes[0].amlList);
        label1 = " \t";
        label2 = " \t";
        #print(genes[0].allList); 
            
        j = 2; 
            
        for i in range(0, expCount):
            
            label1 += "EXP" + str(i + 1) + "\t";
            label2 += "(" + types[j]  + ")\t";
            j += 2;
            
        file.write(label1 + "\n");
        file.write(label2 + "\n");
            
        WriteGenes(genes, file);
        
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nFirst affymetrics save successful.\n");
    
    file.close();
    return 0;


# Save Testing Data, Part II 4.b
def SaveAffymetrics2(genes, types):
    
    fileName = "affymetrics2.txt";
        
    try:
        
        file = open(fileName, "w");
        
        # ----==== Change label second here to however needed. ====----
            
        #expCount = len(genes[0].allList) + len(genes[0].amlList);
        #label = " \t";
            
        #for i in range(0, expCount):
                
            #label += "Test_Exp" + str(i + 1) + "\t";
            
        #file.write(label + "\n");
            
        WriteGenes(genes, file);
        
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nSecond affymetrics save successful.\n");
    
    file.close();
    
    return 0; 

   
# Save Testing Data, Part II 5 
def SaveAffymetrics3(genes, types):
    
    fileName = "affymetrics3.txt";
        
    try:
        
        file = open(fileName, "w");
            
        expCount = len(genes[0].allList) + len(genes[0].amlList);
        label = " \t";
            
        for i in range(0, expCount):
                
            label += "Test_Exp" + str(i + 1) + "\t";
            
        file.write(label + "\n");
            
        WriteGenes(genes, file);
        
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nThird affymetrics save successful.\n");
            
    file.close();
    return 0;



def WriteGenes(genes, file):
    
    for i in range(0, len(genes)):
                
        file.write(genes[i].accessionId);
                
        for k in range(len(genes[i].allList)):
                    
            file.write("\t" + str(genes[i].allList[k][0]));
                
        for k in range(len(genes[i].amlList)):
                    
            file.write("\t" + str(genes[i].amlList[k][0]));
                
        file.write("\n");
    
    return 0;



# Part I, 2.
def Preprocess(genes):
    
    # Eliminate the genes with all As across the experiments, and replace all the expression values below some 
    # threshold cut-off value to that threshold value (pick 20 to be the threshold cut-off value).
    geneCount = len(genes);
    i = 0;
    
    # All genes have the same amount so just set here. 
    allCount = len(genes[0].allList);
    amlCount = len(genes[0].amlList);
    
    while (i < 60): #--------------------------------------------- Set back to geneCount, lowered for faster testing purposes. 
        
        aCount = 0;
        
        for j in range(0, allCount):
            
            if ("A" in genes[i].allList[j][1]):
                ++aCount;
            
        
        for j in range(0, amlCount):
            
            if ("A" in genes[i].amlList[j][1]):
                ++aCount;
        
        # [min, max]
        minMax = GetMinMax(genes[i].allList + genes[i].amlList);
        
        # Eliminate genes with either all A tags, and eliminate the genes with less than two
        # fold change across the experiments, where: max(exp1...exp38)/min(exp1...exp38) < 2. 
        if ( (aCount == allCount + amlCount) or (minMax[1]/minMax[0] < 2) ):
            del genes[i]; 
        else:
            i += 1; 
    
    return 0;


def GetMinMax(list): 
    
    max = 0;
    min = list[0][0];
    
    for i in range(0, len(list)):
        if (list[i][0] > max ):
            max = list[i][0];

        if (list[i][0] < min):
            min = list[i][0];
    
    return [min, max];



# Part II, 4.
def ProcessTraining(genes):
    
    
    return 0;


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
            
        for i in range(0, len(labelsTesting)):
                
            file.write(labelsTesting + "\t");
        
        file.write("\n");
        
        for i in range(0, len(testingPrediction)):
                
            file.write(testingPrediction + "\t");
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nClassification result save successful.\n");
            
    file.close();
    return 0;


def SelectTestingGenes(genesTraining, genesTesting): 
    
    selectedGenes = [];
    
    # For each of the top 50 selected genes in the training data, get the matching ones
    # from the testing data for future classification. 
    for i in range(0, len(genesTraining)):
        
        selectedGenes.append(genesTesting[genesTraining[i].ID]);
    
    return selectedGenes;


def main():
    
    inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
    genesTraining = inputTraining[0];
    typesTraining = inputTraining[1];
    #print(types);
    
    Preprocess(genesTraining); 
    print("\nProcessed.\n");
    SaveAffymetrics1(genesTraining, typesTraining);
    
    # TODO - T test and excel stuff selection of top 50 genes based on p-value
    
    inputTesting = readInData("Loading Leuk_ALL_AML.test .....", "Leuk_ALL_AML.test.txt", 20);
    genesTesting = inputTesting[0];
    typesTesting = inputTesting[1];
    
    # TODO Makes sure only top 50 selected training genes are sent to this function. 
    # Take from the testing data, the matching genes to the selected top 50 training data genes. 
    # genesTesting = SelectTestingGenes(genesTraining, genesTesting);
    
    SaveAffymetrics3(genesTraining, typesTesting);
    
    # ClassifyGenes(genesTraining, typesTraining, genesTesting, typesTesting)
    
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
