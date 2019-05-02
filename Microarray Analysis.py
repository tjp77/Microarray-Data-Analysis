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
                    
                    genes.append(Gene(accession, all, aml, id, allCategory, amlCategory)); #print(genes[-1].allList, "\n\n");
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
                    #print (len(aml), " - ", len(all));
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
        #print(genes[0].allList); 
            
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
    
    #orig = [ [1,2,3], [4,5,6]];
    #print (orig);
    new = [];
    
    i = 0; j = 0; k = 1;
    
    while (i<len(genes[0])):
        
        new.append([genes[0][i]]);
        
        while (k<len(genes)):
            
            new[i].append(genes[k][i]);
            k+=1;
            
        i+=1; k = 1;
    
    #print (new);
    return new;

# Part III, KNN
# https://scikit-learn.org/stable/modules/generated/sklearn.svm.libsvm.fit.html#sklearn.svm.libsvm.fit
def ClassifyGenes(genesTraining, labelsTraining, genesTesting, labelsTesting):
    
    knn = KNeighborsClassifier(algorithm = 'auto', leaf_size = 30,  
                               metric = 'minkowski', metric_params = None, n_jobs = 1,
                               n_neighbors = 3, p = 2, weights = 'uniform');
    #print("classify func genes exp train: ", len(genesTraining)) OK
    #print("classify func genes train: ", len(genesTraining[0])) OK
    #print (labelsTraining) OK
    #print (genesTraining[0]) SEEMS OK
    knn.fit(genesTraining, labelsTraining);  
    
    prediction = knn.predict(genesTesting); 
    
    print("Predic: ", end='');
    #print(prediction);
    for i in range(0, len(prediction)):
        print(prediction[i], end='');
    
    print("\nTarget: ", end='');
    
    for i in range(0, len(labelsTesting)):
        print(labelsTesting[i], end='');
    
    #WriteClassifyResultsToFile(prediction, labelsTesting);
    
    
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


# Select matching testing genes corresponding to traing selections, and form
# label and gene expression arrays in the format required by knn.fit(). 
# make sure only training genes list with top 50, or however many selected genes in passed in, not full genes list.
def SelectTestingGenes(top50, genesTraining, genesTesting, typesTraining, typesTesting): 
    
    selectedTraining = [];
    selectedTesting = [];
    typeToKNNLabel = { "ALL" : 0, "AML" : 1};
    labelsTraining = typesTraining; #len(typesTraining)
    labelsTesting = typesTesting[2:-1:2];
    kNNLabelsTraining = [];
    kNNLabelsTesting = [];
    
    #print (len(labelsTraining))
    #print (len(labelsTesting))
    #print (len(genesTraining[0].allList), " - ", len(genesTraining[0].amlList))
    #print (len(genesTesting[0].allList + genesTesting[0].amlList))
    
    # For each of the top 50 selected genes in the training data, get the matching ones
    # from the testing data for future classification.  
    for i in range(0, len(top50)):
        
        for k in range(0, len(genesTraining)):
            
            if (genesTraining[k].accessionId == top50[i].accessionId):
                selectedTraining.append( np.append(genesTraining[k].allList, genesTraining[k].amlList, axis = 0) );
                selectedTesting.append( np.append(genesTesting[genesTraining[k].ID].allList, genesTesting[genesTraining[k].ID].amlList, axis = 0) ); 
                break;
    
    for i in range(1, len(labelsTraining)):
        
        kNNLabelsTraining.append(typeToKNNLabel[labelsTraining[i][1:4]]);
    
    for i in range(0, len(labelsTesting)):
        
        kNNLabelsTesting.append(typeToKNNLabel[labelsTesting[i][:3]]);
        
    #print (len(selectedTraining[0])) # ---- Coming out ok. 
    
    return [selectedTraining, kNNLabelsTraining, selectedTesting, kNNLabelsTesting];



def DisplaAccuracy(results, targetResults):
    
    correctCount = 0;
    i = 0;
    while (i < len(targetResults)):
        
        if (results[i] == targetResults[i]):
            ++correctCount;
        ++i;
    
    print ("Predicted Results Accuracy: ", correctCount/length, "%")
    
    return 0;


def main():

    if args.program_part == "pre":
        inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
        genesTraining = inputTraining[0];
        typesTraining = inputTraining[1];
        #print(types);
    
        Preprocess(genesTraining); 
        print("\nProcessed.\n");
        SaveAffymetrics1(genesTraining, typesTraining);
    
    elif args.program_part == "post":
        inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
        genesTraining = inputTraining[0];
        typesTraining = inputTraining[1]; #TODO can take this out ince get top 50 working since will use it instead. 
        #print(types);
    
        Preprocess(genesTraining); 

        top_50_input = readInData2("Reading the top 50 genes...", "Affymetrics_top50.txt", 20)
        top_50_genes = top_50_input[0]
        top_50_types = top_50_input[1]
        
        
        #readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);#
        inputTesting = readInData("Loading Leuk_ALL_AML.test .....", "Leuk_ALL_AML.test.txt", 20);
        genesTesting = inputTesting[0];
        typesTesting = inputTesting[1]; #print (len(top_50_genes[0].allList), " - ", len(top_50_genes[0].amlList))
        
        # TODO Makes sure only top 50 selected training genes are sent to this function. 
        # Take from the testing data, the matching genes to the selected top 50 training data genes. 
        genesKNNArrs = SelectTestingGenes(top_50_genes, genesTraining, genesTesting, top_50_types, typesTesting);
        
        #print ( Reformat(genesKNNArrs[0])) # from top 50 post processing on line above
        #print (Reformat(genesKNNArrs[2])) # testing loaded in normally that works for sure
    
        #SaveAffymetrics3(genesTraining, typesTesting); #TODO doesn't seem right, check. 
    
        ClassifyGenes(Reformat(genesKNNArrs[0]), genesKNNArrs[1], Reformat(genesKNNArrs[2]), genesKNNArrs[3]);
        
        
    else: 
        parser.print_help()
        
main();


# ======= To Do ======= 

# Load back in the top50 genes file so can be used in the next part of the program.


# Try dif amounts of n_neighbors for the KNN classifier on full dataset and compare results. 



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

