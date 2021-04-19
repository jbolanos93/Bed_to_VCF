def reference_genome():
    #Coronavirus_2_isolate_genome.fasta
    #Set script to look in reference_genome folder. 
    from pathlib import Path
    path=Path("reference_genome/")
    file_name=''
    for entry in path.iterdir():
        #On MacOS filters the .DS_Store file. 
        if entry.name.startswith("."):
            continue
        else:
            print("*"*100)
            print("\nReference Genome Used:\n '" + entry.name+ "'\n" )
            file_name=entry.name


            #Open File
            fh=open("reference_genome/"+file_name, 'r')
            fn=fh.readlines()
            corona_genome=''
            for line in fn:
                if line.startswith(">"):
                    continue
                else:
                    corona_genome+=line.strip()
    
    return corona_genome
    
def read_bed_file():
    import pandas as pd
    from pathlib import Path
    file_list=[]
    path=Path("bed/")
    file_name=''

    for entry in path.iterdir():
        
        
        file_name=entry.name
        save_file= file_name.strip(".bed")
        bed_file=open("bed/"+file_name, "r")
        
        #On MacOS filters the .DS_Store file. 
        if file_name.startswith("."):
            continue
        else: 
            
            bed_file_list=bed_file.readlines()

            yield bed_file_list,save_file
        
        
        #print("bed file found :'" + entry.name+"'" )
        #file_list.append(entry.name)

def create_dataframe(a_list,corona_genome):
    
    import pandas as pd
    from pathlib import Path
    root="output_vcf/"
    new_file=[]
    print("*"*100)
    print("Converting File: "+ a_list[1]+".bed\n")
    print("Indels found: ")
    path=Path("bed/") 
    for line in a_list[0]:

        df_list=line.split("\t")
        ID=df_list[len(df_list)-1].strip()
        deletion_start=int(df_list[len(df_list)-3].strip())-1
        deletion_end=int(df_list[len(df_list)-2].strip())
        #print(deletion_start,deletion_end)
        if ID.startswith("d") or ID.startswith("i"):


            #Sets pos equal to pos-1
            #to follow the sample vcf file. 
            b=line.split()
            b[2]=b[1]
            b='\t'.join(b)
            #print(b)


            original=b.strip()


            original+="\t" +corona_genome[deletion_start:deletion_end]+"\t"+corona_genome[deletion_start]

            new_file.append(original)



            print(" •",ID[len(ID)-len(ID):] ,": ",corona_genome[deletion_start:deletion_end])


        else:    
            original=line.strip()
            original+="\t" +ID[len(ID)-len(ID)]+"\t"+ID[len(ID)-1]
            new_file.append(original)

        
    df = pd.DataFrame([sub.split("\t") for sub in new_file],index=None,columns=["#CHROM", "POS-1", "POS", 'ID',"REF", "ALT"])
    df.drop('POS-1', inplace=True, axis=1)
    
    df.to_csv(root+a_list[1]+".vcf",sep="\t",header=True, index=False)    
    print("\nPreview:  " +a_list[1]+ ".vcf\n\n",df,"\n")


def main():
    import pandas as pd
    from pathlib import Path
    
    print("Converting Files: ")
    path=Path("bed/")
    
    for entry in path.iterdir():
        if entry.name.startswith("."):
            continue
        else:
            print(entry.name)
        
    
    reference=reference_genome()
    
    #Use Generator
    for i in read_bed_file():
        create_dataframe(i,reference)
    
    print(".vcf files saved in '/output_vcf' directory.")

if __name__ == "__main__":
    main()

    


    

