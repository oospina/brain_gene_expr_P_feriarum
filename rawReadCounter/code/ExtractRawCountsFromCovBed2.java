import java.io.*;
import java.util.*;

//THIS PROGRAM IDENTIFIES THE NUMBER OF READS MAPPING TO EACH GENE, COLLAPSING ISOFORMS BUT ONLY COUNTING EACH READ ONCE 
//NOTE: NAMES OF READS MAPPED TO MULTIPLE GENES NOT COUNTED

//EXAMPLE USAGE: java -Xmx16g ExtractRawCountsFromCovBed2 I29357_aligned_SORTED_transcript.cov.txt I29357_RawCountsNew
public class ExtractRawCountsFromCovBed2{

	public static void main(String[] args){
		try{

			String inFile = args[0];							//fasta or fastq file
			String outStem = args[1];

			//SETUP HASHMAP FOR COUNTING READS MAPPED ONE ENTRY PER GENE WITH GENE NAME AS KEY
			//HashMap<String,Integer>coverageMap = new HashMap<String,Integer>();	

			//SETUP HASHMAP FOR IDENTIFYING GENE TO WHICH EACH READ IS MAPPED. READS MAPPING TO MULTIPLE GENES WILL GET VALUE MULTIPLE NAD CAN BE DISQUALIFIED DOWNSTREAM
			HashMap<String,String>readMap = new HashMap<String,String>();

			BufferedReader br = new BufferedReader(new FileReader(new File(inFile)));
			String tempS=br.readLine(); //e.g. TRINITY_DN302420_c0_g1_i1	64	186	A00327:34:HG23NDSXX:4:1401:31141:22420	60	+	64	186	0,0,0	1	122,	0,	0	0122	0.0000000
			String gene="";
			String iso="";
			String geneMapped="";
			String read="";
			Integer count=0;

			//int nMappingsToMultipleIsoforms=0;
			//int nMappingsToMultipleGenes=0;
			//int nMappingsTotal=0;
			int nReads=0;
			int nGenes=0;
			int nIsoforms=0;
			int nReadsMappedToSingleGene=0;

			//int nUnmappedReads=0;
			//int nReadsMappedToMultipleGenes=0;

			int nLines=0;

			while(tempS!=null){

				gene=tempS.split("\t")[0].split("_i")[0];

				if(gene.equals(".")){
					tempS=br.readLine();
					continue;	//allow unmapped reads to be dispersed throughout list
					//break;		//assume all unmapped reads are listed at end of file
				}

				nLines++;
				if(nLines%10000==0){System.out.print("\r"+nLines);}

				read=tempS.split("\t")[3];
				geneMapped=readMap.get(read);
				if(geneMapped==null){			//first mapping of this read
					nReads++;
					readMap.put(read,gene);
				}else if(!geneMapped.equals(gene)){
					readMap.put(read,"MULTIPLE");	//read mapped to multiple genes, stop counting
				}

				tempS=br.readLine();
			}
			br.close();
			System.out.println();
			System.out.println(nReads+" reads were observed.");			
			System.out.println(nLines+" mappings were observed.");

			//GET COUNTS FOR READS THAT ONLY MAP TO ONE GENE(BUT ALLOW MULTIPLE ISOFORMS)
			HashMap<String,Integer>geneMap = new HashMap<String,Integer>();	
			HashMap<String,Integer>isoformMap = new HashMap<String,Integer>();	

			br = new BufferedReader(new FileReader(new File(inFile)));
			tempS=br.readLine();
			nLines=0;
			while(tempS!=null){
				nLines++;
				if(nLines%10000==0){System.out.print("\r"+nLines);}

				gene=tempS.split("\t")[0].split("_i")[0];
				iso=tempS.split("\t")[0];
				read=tempS.split("\t")[3];


				if(gene.equals(".") || gene.equals("*")){
					break;		//assume all unmapped reads are listed at end of file
				}

				if(readMap.get(read).equals("MULTIPLE")){	//multiple genes, do not count at all
					tempS=br.readLine();
					continue;
				}


				geneMapped=readMap.get(read);

				//map to isoform level
				count=(Integer)isoformMap.get(iso);
				if(count==null){
					nIsoforms++;
					count=0;					
				}
				isoformMap.put(iso,(Integer)(count+1));

				if(readMap.get(read).equals("MULTIPLE2")){	//already counted for this gene, don't count again
					tempS=br.readLine();
					continue;
				}

				nReadsMappedToSingleGene++;

				count=(Integer)geneMap.get(gene);
				if(count==null){
					nGenes++;
					count=0;
				}
				geneMap.put(gene,(Integer)(count+1));

				readMap.put(read,"MULTIPLE2");	//disable read so that it does not get counted again for this gene

				tempS=br.readLine();
			}
			br.close();
			System.out.println();
			System.out.println(nGenes+" genes had at least one good read mapped.");
			System.out.println(nIsoforms+" isoforms had at least one good read mapped.");
			System.out.println(nReadsMappedToSingleGene+" reads were mapped to a single gene.");

			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outStem+"_genes.csv")));
			Iterator it = geneMap.entrySet().iterator();
			while (it.hasNext()) {
				Map.Entry entry = (Map.Entry) it.next();
				gene=(String)entry.getKey();
				count= (Integer)entry.getValue();
				bw.write(gene+","+count+"\n");
			}
			bw.flush();
			bw.close();


			bw = new BufferedWriter(new FileWriter(new File(outStem+"_isoforms.csv")));
			it = isoformMap.entrySet().iterator();
			while (it.hasNext()) {
				Map.Entry entry = (Map.Entry) it.next();
				iso=(String)entry.getKey();
				count= (Integer)entry.getValue();
				bw.write(iso+","+count+"\n");
			}
			bw.flush();
			bw.close();

		}catch(IOException ioe){System.out.println(ioe.getMessage());}
	}

}
