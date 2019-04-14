#!/user/bin/python


def readInData(prompt):
    
    data = "";
    not_read = True
    
    while (not_read):
        
        try:
            
            file_name = input(prompt);
            
            with open(file_name, "r") as file:
                for line in file:
                    
                    # strip to be safe due to proir input whitespace problems.     
                    data += line.strip();
        except:
        
            print("\nTrouble reading input file. Make sure name is correct and the file is in the right format.\n");
            
        else:
            print("\nFile read successfully.\n");
            not_read = False;
            file.close();
    
    return data.strip();



# Part I, 2.
def Preprocess()
    
    
    return 0;


# Part II, 5.
def Process()
    
    
    return 0;


def main()
    
    # plh
    
    return 0;
    


main();


# ======= Notes =======
#
# "A class discovery procedure automatically discovered the distinction between
# acute myeloid leukemia (AML) and acute lymphoblastic leukemia (ALL)"
#
#
