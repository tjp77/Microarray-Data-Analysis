#!/user/bin/python


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
        
        try:
            
            fileName = input(prompt); 
            
            with open(fileName, "r") as file:
                
                # strip any input to be safe due to proir experience with input whitespace problem.
                firstLine = file.readline().split("\t").strip(); 
                
                # Not sure will need these, but have to pass through them either way. 
                line = file.readline(); 
                line = file.readline(); 
                
                
                # Go through remaining lines and get all things interested in.
                while line:
                    
                    line = file.readline().split("\t").strip(); 
                    
                     # For all where first col, description contains "(endogenous control)", skip since supossed to eliminate them.
                    if ("endogenous control" in line[0]):
                        continue;
                    
                    # second column is Accession ID.
                    gene - Gene(line[1]);
                    
                    j = 2;
                    while (j < len(line)):
                        
                        if ("ALL" in firstLine[j]):
                            
                            # ex: [88, A] for line 7 of input data file
                            gene.allList.append([line[j], line[j+1]]);
                            
                        elif ("AML" in firstLine[j]):
                            
                            gene.amlList.append([line[j], line[j+1]]);
                        
                        j += 2;
                    
                    
                    genes.append(gene);

        except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            
        else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return genes;



def SaveAffymetrics1():
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetrics1.txt";
        
        try:
        
            file = open(fileName, "w");
            
            #file.write(...);
        
        except:
        
            print("\nError writing to file.\n");

        else:
            print("\nTable save successful.\n");
            notSaved = False;
    
    file.close();
    return 0;
    
    
def SaveAffymetrics2():
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetric2s.txt";
        
        try:
        
            file = open(fileName, "w");
            
            #file.write(...);
        
        except:
        
            print("\nError writing to file.\n");

        else:
            print("\nTable save successful.\n");
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
def Process():
    
    
    return 0;



# Par III, KNN
def classify():
    
    
    return 0;


def main():
    
    genes = readInData("Enter then name of the file to load:\n");
    
    return 0;
    


main();

# ======= To Do =======

# - Output read in file and print to test if coming in ok. 


# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#
