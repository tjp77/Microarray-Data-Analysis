#!/user/bin/python
import sys;
import xlsxwriter
import numpy as np
from sklearn import datasets
from sklearn.neighbors import KNeighborsClassifier

workbook = xlsxwriter.Workbook('microarray-analysis.xlsx')
affy1_worksheet = workbook.add_worksheet('affymetrics1')

class Gene:
    
    def __init__(self, acc, all, aml):
        
        self.accessionId = acc;
        self.allList = all;
        self.amlList = aml;


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
                    while (j < len(line) - 1):
                        
                        if ( (int)(line[j]) < threshold ):
                            line[j] = threshold; 
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            all.append([(int)(line[j]), line[j+1]]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            aml.append([(int)(line[j]), line[j+1]]);
                        
                        j += 2;
                    
                    genes.append(Gene(accession, all, aml)); #print(genes[-1].allList, "\n\n");
                    line = file.readline(); 
        except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");p
            sys.exit();
            
        else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return [genes, firstLine];


# Save Training Data, Part I
def SaveAffymetrics1(genes, types):
    expCount = len(genes[0].allList) + len(genes[0].amlList);
    label1 = ""
    label2 = "";

    j = 2

    row = 0
    col = 1

    for i in range(expCount):
        label1 = "EXP" + str(i + 1)
        label2 = "(" + types[j]  + ")"
        j += 2;

        affy1_worksheet.write(row, col, label1)
        affy1_worksheet.write(row+1, col, label2)

        col += 1

    WriteGenes(genes, affy1_worksheet)
    
    # fileName = "affymetrics1.txt";
        
    # try:
            
    #     file = open(fileName, "w");
            
    #     expCount = len(genes[0].allList) + len(genes[0].amlList);
    #     label1 = ""
    #     label2 = "";
    #     #print(genes[0].allList); 
            
    #     j = 2; 

    #     row = 0
    #     col = 1
            
    #     for i in range(0, expCount):
            
    #         label1 = "EXP" + str(i + 1)
    #         label2 = "(" + types[j]  + ")"
    #         j += 2;
    #         affy1_worksheet.write(row, col, label1)
    #         affy1_worksheet.write(row+1, col, label2)
    #         col += 1
            
    #     file.write(label1 + "\n");
    #     file.write(label2 + "\n");

            
    #     WriteGenes(genes, file);
        
    # except:
        
    #     print("\nError writing to ", fileName, "file.\n");

    # else:
    #     print("\nFirst affymetrics save successful.\n");
    
    # file.close();
    # return 0;


# Save Testing Data, Part II 4.b
def SaveAffymetrics3(genes, types):
    
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



def WriteGenes(genes, worksheet):

    row = 2
    col = 0
    for i in range(0, len(genes)):
        col = 0
        worksheet.write(row, col, genes[i].accessionId)
        col += 1
                
        for k in range(len(genes[i].allList)):
            worksheet.write(row, col, str(genes[i].allList[k][0]))
            col += 1

        for k in range(len(genes[i].amlList)):
            worksheet.write(row, col, str(genes[i].amlList[k][0]))
            col += 1
                
        row += 1
    
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

    # FIXME: Set upper limit to size of dataset
    while (i < 5): #--------------------------------------------- Set back to geneCount, lowered for faster testing purposes. 
        
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


# Par III, KNN
def classifygenes(genes):
    
    knnSize = 0.6 * len(genes);
    
    knn = KNeighborsClassifier(algorithm = 'auto', leaf_size = 30, # what was leaf size? 
                               metric = 'minkowski', metric_params = None, n_jobs = 1,
                               n_neighbors = 3, p = 2, weights = 'uniform');
    
    trainData = "";
    trainTargetLabels = "";
    
    testData = "";
    testTargetLabels = "";
    
    knn.fit(trainData, trainTargetLabels);
    
    print("Classifier predictions for testing data:");
    print(knn.predict(testData));
    
    print("Target labels for testing data:");
    print(testTargetLabels);
    
    return 0;


def main():
    
    # Part I
    inputTraining = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....", "ALL_vs_AML_train_set_38_sorted.txt", 20);
    genesTraining = inputTraining[0];
    typesTraining = inputTraining[1];
    #print(types);
    
    Preprocess(genesTraining); 
    print("\nProcessed.\n");
    SaveAffymetrics1(genesTraining, typesTraining);
    
    
    # Part II
    inputTesting = readInData("Loading Leuk_ALL_AML.test .....", "Leuk_ALL_AML.test.txt", 20);
    genesTesting = inputTesting[0];
    typesTesting = inputTesting[1];
    
    SaveAffymetrics3(genesTesting, typesTesting);
    
    
    return 0;
    


main();

workbook.close()

# ======= To Do ======= 

# PreProcess function:
# Eliminate the genes with less than two fold change across the experiments (max/min <2);
# - REVIEW SAVED AFFY FILE, SEE IF SEEMS LIKE RIGHT THINGS REMOVED. ('seems' due to file too large to check all)



# Sort genes by p-values/T test stuff Part II 4. a


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
