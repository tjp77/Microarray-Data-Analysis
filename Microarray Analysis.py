#!/user/bin/python
import sys;
import numpy as np
import argparse
from sklearn import datasets
from sklearn.neighbors import KNeighborsClassifier

parser = argparse.ArgumentParser(description="Run code for the Microarray project.")
parser.add_argument("program_part", metavar="pre | post", type=str, help="To run pre-Excel functions use pre, to run post-excel functions use post")
args = parser.parse_args()

class Gene:
    
    def __init__(self, acc, all, aml, rowID, allCategory, amlCategory):
        
        self.accessionId = acc;
        self.allList = all;
        self.amlList = aml;
        self.allCategoryList = allCategory;
        self.amlCategoryList = amlCategory;
        self.ID = rowID;


def readInData(prompt, fileName, threshold=20):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        try:
            
            print(prompt); 
            
            with open(fileName, "r") as file:
                
                firstLine = file.readline().split("\t"); 
                
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
                    allCategory = [];
                    amlCategory = [];
                    
                    j = 2;
                    # Get each expresson value.
                    while (j < len(line) - 1):
                        
                        if ( (int)(line[j]) < threshold ):
                            line[j] = threshold; 
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            all.append((int)(line[j]));
                            allCategory.append(line[j+1]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            aml.append((int)(line[j]));
                            amlCategory.append(line[j+1]);
                        
                        j += 2;
                    
                    genes.append(Gene(accession, all, aml, id, allCategory, amlCategory)); 
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


def readInData2(prompt, fileName, threshold=20):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        try:
            
            print(prompt); 
            
            with open(fileName, "r") as file:
                
                file.readline();
                firstLine = file.readline().split("\t"); 
                
                # Keep track of which file line gene is on, so can easily pull testing genes corresponding to selected top 50 training genes. 
                id = 0;
                line = file.readline();
                
                # Go through remaining lines and get all things interested in.
                while line:
                    
                    line = line.split("\t")
                        
                    # second column is Accession ID.
                    accession = line[0];
                    all = [];
                    aml = [];
                    
                    j = 1;
                    # Get each expresson value.
                    while (j < len(firstLine)):
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            all.append((int)(line[j]));
                            
                        elif ("AML" in firstLine[j]):
                            
                            aml.append((int)(line[j]));
                        
                        j += 1;
                    
                    genes.append(Gene(accession, all, aml, id, [], [])); 
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
            
        j = 2; 
            
        for i in range(0, expCount):
            
            label1 += "EXP" + str(i + 1) + "\t";
            label2 += "(" + types[j][:3]  + ")\t";
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
                    
            file.write("\t" + str(genes[i].allList[k]));
                
        for k in range(len(genes[i].amlList)):
                    
            file.write("\t" + str(genes[i].amlList[k]));
                
        file.write("\n");
    
    return 0;



# Part I, 2.
def Preprocess(genes):
    
    # Eliminate the genes with all As across the experiments, and replace all the expression values below some 
    # threshold cut-off value to that threshold value (pick 20 to be the threshold cut-off value).
    geneCount = len(genes);
    i = 0;
    
    # All genes have the same amounts so just set here. 
    allCount = len(genes[0].allList);
    amlCount = len(genes[0].amlList);
    
    while (i < len(genes)):  
        
        aCount = 0;
        
        for j in range(0, allCount):
            
            if ("A" in genes[i].allCategoryList[j]):
                ++aCount;
            
        for j in range(0, amlCount):
            
            if ("A" in genes[i].amlCategoryList[j]):
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
    min = list[0];
    
    for i in range(0, len(list)):
        if (list[i] > max ):
            max = list[i];

        if (list[i] < min):
            min = list[i];
    
    return [min, max];

# Reformat genes array from that it was loaded in as to the format knn expects. 
def Reformat(genes):
    
    new = [];
    
    i = 0; j = 0; k = 1;
    
    while (i<len(genes[0])):
        
        new.append([genes[0][i]]);
        
        while (k<len(genes)):
            
            new[i].append(genes[k][i]);
            k+=1;
            
        i+=1; k = 1;

    return new;

# Part III, KNN
# https://scikit-learn.org/stable/modules/generated/sklearn.svm.libsvm.fit.html#sklearn.svm.libsvm.fit
def ClassifyGenes(genesTraining, labelsTraining, genesTesting, labelsTesting):
    
    k = 5;
    alg = "auto";
    metric = "minkowski";
    weight = "uniform";
    knn = KNeighborsClassifier(algorithm = alg, leaf_size = 30,  
                               metric = metric, metric_params = None, n_jobs = 1,
                               n_neighbors = k, p = 2, weights = weight);
    
    knn.fit(genesTraining, labelsTraining);  
    prediction = knn.predict(genesTesting); 
    
    total = len(prediction);
    correctMatchCount = 0;
    predictedResults = "Predic: ";
    targetResults = "Target: ";
    
    for i in range(0, len(prediction)):
        
        predictedResults += str(prediction[i]); 
        targetResults += str(labelsTesting[i]);
        
        if (str(prediction[i]) == str(labelsTesting[i])):
            correctMatchCount += 1;
    
    print (predictedResults);
    print (targetResults);
    
    print ("\nPredicted Results Accuracy: ", str(round(correctMatchCount/total, 3) * 100) + "%\nWith:");
    print ("k = " + str(k), "\nMAlgorithm = " + alg, "\nMetric = " + metric, "\n" + weight + " weight");
    
    WriteClassifyResultsToFile(prediction, labelsTesting);
    
    return 0;


def WriteClassifyResultsToFile(testingPrediction, labelsTesting):
    
    fileName = "classificationResults.txt";
        
    try:
        
        file = open(fileName, "w");
            
        for i in range(0, len(labelsTesting)):
                
            file.write(str(labelsTesting[i]) + "\t");
        
        file.write("\n");
        
        for i in range(0, len(testingPrediction)):
                
            file.write(str(testingPrediction[i]) + "\t");
    except:
        
        print("\nError writing to ", fileName, "file.\n");

    else:
        print("\nClassification result save successful.\n");
        
    file.close();
    return 0;


# Select matching testing genes corresponding to traing selections, and form
# label and gene expression arrays in the format required by knn.fit(). 
# make sure only training genes list with top 50, or however many selected genes in passed in, not full genes list.
def SelectTestingGenes(top50, genesTesting, top50Types, typesTesting): 
    
    selectedTraining = [];
    selectedTesting = [];
    typeToKNNLabel = { "ALL" : 0, "AML" : 1};
    labelsTraining = top50Types; 
    labelsTesting = typesTesting[2:-1:2];
    kNNLabelsTraining = [];
    kNNLabelsTesting = [];
    
    # For each of the top 50 selected genes in the training data, get the matching ones
    # from the testing data for future classification.  
    for i in range(0, len(top50)):
        
        for k in range(0, len(genesTesting)):
            
            if (genesTesting[k].accessionId == top50[i].accessionId):
                
                selectedTraining.append( np.append(top50[i].allList, top50[i].amlList, axis = 0) );
                selectedTesting.append( np.append(genesTesting[k].allList, genesTesting[k].amlList, axis = 0) ); 
                break;
    
    for i in range(1, len(labelsTraining)):
        
        kNNLabelsTraining.append(typeToKNNLabel[labelsTraining[i][1:4]]);
    
    for i in range(0, len(labelsTesting)):
        
        kNNLabelsTesting.append(typeToKNNLabel[labelsTesting[i][:3]]);
        
    return [selectedTraining, kNNLabelsTraining, selectedTesting, kNNLabelsTesting];


def main():

    if args.program_part == "pre":
        
        inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
        genesTraining = inputTraining[0];
        typesTraining = inputTraining[1];
    
        Preprocess(genesTraining); 
        print("\nProcessed.\n");
        SaveAffymetrics1(genesTraining, typesTraining);
        
        # SaveAffymetrics2 would have been the second file save requested in instructions, the
        # post T-Test file, had it been done in code rather then excel. Numbers kept as if 
        # all three were done here for easier tracking of what was from what step. 
    
    elif args.program_part == "post":

        top_50_input = readInData2("Reading the top 50 genes...", "tmp_top.txt", 20)
        top_50_genes = top_50_input[0]
        top_50_types = top_50_input[1]
        
        inputTesting = readInData("Loading Leuk_ALL_AML.test .....", "Leuk_ALL_AML.test.txt", 20);
        genesTesting = inputTesting[0];
        typesTesting = inputTesting[1]; 
        
        SaveAffymetrics3(genesTesting, typesTesting);
        
        # Make sure only top 50 selected training genes are sent to this function. 
        # Take from the testing data, the matching genes to the selected top 50 training data genes. 
        genesKNNArrs = SelectTestingGenes(top_50_genes, genesTesting, top_50_types, typesTesting);
        
        ClassifyGenes(Reformat(genesKNNArrs[0]), genesKNNArrs[1], Reformat(genesKNNArrs[2]), genesKNNArrs[3]);
        
    else: 
        parser.print_help()
        
main();


# ======= To Do ======= 


# Try dif amounts of n_neighbors for the KNN classifier on full dataset and compare results. 



# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#

