#!/user/bin/python


def readInData(prompt):
    
    data = "";
    notRead = True
    
    while (notRead):
        
        try:
            
            fileName = input(prompt);
            
            with open(fileName, "r") as file:
                for line in file:
                    
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
