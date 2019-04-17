#!/user/bin/python
import sys;

class Gene:
    
    allList = [];
    amlList = [];
    accessionId = "";
    
    def __init__(acc):
        
        accession = acc;


def readInData(prompt):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        #try:
            
            fileName = "ALL_vs_AML_train_set_38_sorted.txt";
            print(prompt); 
            
            with open(fileName, "r") as file:
                
                firstLine = file.readline().split("\t"); 
                
                #for i in range(0, len(firstLine)):
                 #   print (str(i) + firstLine[i] + "|");
                
                # Not sure will need these, but have to pass through them either way. 
                line = file.readline(); 
                line = file.readline(); 
                
                
                # Go through remaining lines and get all things interested in.
                while line:
                    
                    line = file.readline().split("\t")
                     # For all where first col, description contains "(endogenous control)", skip since supossed to eliminate them.
                    if ("endogenous control" in line[0]):
                        continue;
                    
                    # second column is Accession ID.
                    gene = Gene(line[1]);
                    
                    j = 2;
                    while (j < len(line)):
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            gene.allList.append([line[j], line[j+1]]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            gene.amlList.append([line[j], line[j+1]]);
                        
                        j += 2;
                    
                    
                    genes.append(gene);

        #except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            
        #else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return [genes, firstLine];



def SaveAffymetrics1(genes, types):
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetrics1.txt";
        
        try:
        
            file = open(fileName, "w");
            
            expCount = len(genes[0].allList);
            label1 = " \t";
            label2 = " \t";
            
            j = 2;
            for i in range(1, expCount):
                
                label1 += "EXP" + i + "\t";
                label2 += "(" + types[j]  + ")\t";
                j += 2;
            
            file.write(label1 + "\n");
            file.write(label2 + "\n");
            
            for i in range(0, genes):
                
                file.write(genes[i].accession);
                
                for k in range(len(genes[i].allList)):
                    
                    file.write("\t" + genes[i].allList[k][0]);
                
                for k in range(len(genes[i].amlList)):
                    
                    file.write("\t" + genes[i].amlList[k][0]);
        
        except:
        
            print("\nError writing to file.\n");

        else:
            print("\nTable save successful.\n");
            notSaved = False;
    
    file.close();
    return 0;
    
    
def SaveAffymetrics2(genes, types):
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetrics2.txt";
        
        try:
        
            file = open(fileName, "w");
            
            file.write(...);
        
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
    
    while (i < geneCount):
        
        aCount = 0;
        allCount = len(genes[i].allList);
        amlCount = len(genes[i].amlList);
        
        for j in range(0, allCount):
            
            if ("A" in genes[i].allList[j][1]):
                ++aCount;
            
            if (genes[i].allList[j][0] < 20):
                genes[i].allList[j][0] = 20;
        
        for j in range(0, amlCount):
            
            if ("A" in genes[i].amlList[j][1]):
                ++aCount;
                
            if (genes[i].amlList[j][0] < 20):
                genes[i].amlList[j][0] = 20;
        
        if (aCount == allCount + amlCount):
            del genes[i];
        else:
            ++i;
    

    # Eliminate the genes with less than two fold change across the experiments (max/min <2);
    
    return 0;


# Part II, 5.
def Process(genes):
    
    
    return 0;



# Par III, KNN
def classifygenes(genes):
    
    
    return 0;


def main():
    
    input = readInData("Enter then name of the file to load:\n");
    genes = input[0];
    types = input[1];
    print(types);
    
    # Part I
    Preprocess(genes)
    SaveAffymetrics1(genes, types);
    
    
    # Part II
    
    
    SaveAffymetrics2(genes, types);
    
    
    return 0;
    


main();

# ======= To Do =======

# - Output read in file and print to test if coming in ok. 
#   - Types being filled ok. 

# - Test/review SaveAffymetrics1 produced data compared to project example out put shown. 

# --- PreProcess function:
# Eliminate the genes with less than two fold change across the experiments (max/min <2);


# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#
