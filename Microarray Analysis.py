#!/user/bin/python
import sys;

class Gene:

    
    def __init__(self, acc, all, aml):
        
        self.accessionId = acc;
        self.allList = all;
        self.amlList = aml;


def readInData(prompt):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        try:
            
            fileName = "ALL_vs_AML_train_set_38_sorted.txt";
            print(prompt); 
            
            with open(fileName, "r") as file:
                
                firstLine = file.readline().split("\t"); 
                
                #for i in range(0, len(firstLine)):
                 #   print (str(i) + firstLine[i] + "|");
                
                # Not sure will need these, but have to pass through them either way. 
                line = file.readline(); 
                line = file.readline(); 
                
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
                    while (j < len(line)):
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            all.append([(int)(line[j]), line[j+1]]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            aml.append([(int)(line[j]), line[j+1]]);
                        
                        j += 2;
                    
                    genes.append(Gene(accession, all, aml)); #print(genes[-1].allList, "\n\n");
                    line = file.readline(); 

        except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            sys.exit();
            
        else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return [genes, firstLine];


# Save Training Data Set
def SaveAffymetrics1(genes, types):
    
    notSaved = True;
    
    while (notSaved):
    
            fileName = "affymetrics1.txt";
        
        #try:
        
            file = open(fileName, "w");
            
            expCount = len(genes[0].allList) + len(genes[0].amlList);
            label1 = " \t";
            label2 = " \t";
            #print(genes[0].allList);
            
            j = 2; 
            
            for i in range(1, expCount):
                #print(j, " - ", i, "\n");
                label1 += "EXP" + str(i) + "\t";
                label2 += "(" + types[j]  + ")\t";
                j += 2;
            
            file.write(label1 + "\n");
            file.write(label2 + "\n");
            
            for i in range(0, len(genes)):
                
                file.write(genes[i].accessionId);
                
                for k in range(len(genes[i].allList)):
                    
                    file.write("\t" + str(genes[i].allList[k][0]));
                
                for k in range(len(genes[i].amlList)):
                    
                    file.write("\t" + str(genes[i].amlList[k][0]));
                
                file.write("\n");
        
        #except:
        
            print("\nError writing to file.\n");

       # else:
            print("\nTable save successful.\n");
            notSaved = False;
    
    file.close();
    return 0;
    
    
# Save Testing Data Set
def SaveAffymetrics2(genes, types):
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetrics2.txt";
        
        try:
        
            file = open(fileName, "w");
            
            #file.write(...);
        
        except:
        
            print("\nError writing to file.\n");

        else:
            print("\nAffymetrics save successful.\n");
            notSaved = False;
    
    file.close();
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
        
        print(i, "\n");
        aCount = 0;
        
        for j in range(0, allCount):
            
            if ("A" in genes[i].allList[j][1]):
                ++aCount;
            
            # print(genes[i].allList[j][0]);
            if (genes[i].allList[j][0] < 20):
                genes[i].allList[j][0] = 20;
        
        for j in range(0, amlCount):
            
            if ("A" in genes[i].amlList[j][1]):
                ++aCount;
            
            if (genes[i].amlList[j][0] < 20):
                genes[i].amlList[j][0] = 20;
        
        # [min, max]
        minMax = GetMinMax(genes[i].allList + genes[i].amlList);
        
        # Eliminate genes with either all A tags, and eliminate the genes with less than two
        # fold change across the experiments, where: max(exp1...exp38)/min(exp1...exp38) < 2. 
        if ( (aCount == allCount + amlCount) or (minMax[1]/minMax[0] < 2) ):
            del genes[i]; print ("deleted ", i);
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



# Part II, 5.
def Process(genes):
    
    
    return 0;



# Par III, KNN
def classifygenes(genes):
    
    
    return 0;


def main():
    
    input = readInData("Loading ALL_vs_AML_train_set_38_sorted.txt .....");
    genes = input[0];
    types = input[1];
    #print(types);
    
    # Part I
    Preprocess(genes); 
    print("\nProcessed.\n");
    SaveAffymetrics1(genes, types);
    print("\nFirst affymetrics Saved..\n");
    
    
    # Part II
    
    
    #SaveAffymetrics2(genes, types);
    
    
    return 0;
    


main();

# ======= To Do ======= 

# PreProcess function:
# Eliminate the genes with less than two fold change across the experiments (max/min <2);
# - REVIEW SAVED AFFY FILE, SEE IF SEEMS LIKE RIGHT THINGS REMOVED. ('seems' due to file too large to check all)

# Save testing data set function

# Sort genes by p-values/T test stuff Part II 4. a

# Part II 5. 


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
