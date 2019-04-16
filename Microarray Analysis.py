#!/user/bin/python


class Gene:
    
    allList = [];
    amlList = [];
    accession = "";


def readInData(prompt):
    
    # array of gene object, each of gene class, each one line/row of the file data. 
    genes = [];
    notRead = True
    
    while (notRead):
        
        try:
            
            fileName = input(prompt); 
            
            with open(fileName, "r") as file:
                
                line = file.readline(); 
                
                # parse first line, go through find what index turns to AML from ALL. SHould be 58 here for test data file. 
                
                line = file.readline(); 
                
                # Not sure what this line is yet.
                
                line = file.readline(); 
                
                # Not sure what this line is yet, either.
                
                # Go through remaining lines and get all things interested in.
                while line:
                    
                    # File columns tab delimited.
                    
                    # For all where first col, description contains "(endogenous control)", skip since says to eliminate them.
                    
                    # second column is Accession ID.
                    
                    # strip any input to be safe due to proir input whitespace problems.     
                    data += line.strip();
        except:
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            
        else:
            print("\nFile read successfully.\n");
            notRead = False;
            file.close();
    
    return data.strip();



def SaveAffymetrics():
    
    notSaved = True;
    
    while (notSaved):
    
        fileName = "affymetrics.txt";
        
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
def Preprocess():
    
    # Eliminate the genes with all As across the experiments;
    
    # Replace all the expression values below some threshold cut-off value to that threshold value (pick 20 to be the threshold cut-off value);

    # Eliminate the genes with less than two fold change across the experiments (max/min <2);
    
    return 0;


# Part II, 5.
def Process():
    
    
    return 0;


def main():
    
    # plh
    
    return 0;
    


main();


# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#
