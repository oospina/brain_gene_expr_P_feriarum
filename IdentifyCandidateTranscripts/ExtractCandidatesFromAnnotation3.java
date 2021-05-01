import java.io.*;
import java.util.*;

public class ExtractCandidatesFromAnnotation3{
  public static void main(String[] args){
      try{

		String annotationFile = args[0];				//e.g. ref_transcriptome_annotations.txt
		String keywordsFile = args[1];					//e.g. CandidateGeneKeywords.txt
		String outStem = args[2];						//e.g. Candidate
		int nExpected=10000;

		//determine nubmer of keywords
		int nTargets=0;
		String tempS="";
		BufferedReader br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(keywordsFile) ) ));
		tempS=br.readLine();
		while(tempS!=null){
			nTargets++;
			tempS=br.readLine();
		}
		br.close();

		String targets[] = new String[nTargets];
		String keywords[][] = new String[nTargets][50];		//up to 10 keywords allowed per target
		int type[] = new int[nTargets];						//0=>AND, 1=>OR
		int nKeywords[] = new int[nTargets];

		br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(keywordsFile) ) ));
		for(int i=0; i<nTargets; i++){
			tempS=br.readLine();
			targets[i]=tempS;
			tempS=tempS.toUpperCase();
			String tempA[]=tempS.split("\t");	//e.g. AND\t=Syntaxin-\t!=Syntaxin-binding
			if(tempA[0].equals("AND")){type[i]=0;}else if(tempA[0].equals("OR")){type[i]=1;}else{System.out.println("UNRECOGNIZED TYPE "+tempA[0]); return;}
			System.out.print(type[i]);
			for(int j=1; j<tempA.length; j++){
				nKeywords[i]++;
				System.out.print("[");
				keywords[i][j-1]=tempA[j].toUpperCase();
				System.out.print(keywords[i][j-1]);
				System.out.print("]");
			}
			System.out.println();
		}

		String candidateNames[]=new String[nExpected];
		String candidateAnnot[] = new String[nExpected];
		String candidateTargetsFound[] = new String[nExpected];
		int candidateKeywordIDs[] = new int[nExpected];
		int candidateKeywordPos[] = new int[nExpected];
		int candidateCounts[][] = new int[nExpected][nTargets];
		int targetCounts[] = new int[nTargets];		
		int targetCountsNew[] = new int[nTargets];		


		int n;
		String TEMPS="";
		int nCandidates=0;
		int candID=0;
		int nHits=0;
		int nNewHits=0;
		int nHitsTotal=0;
		int nNewHitsTotal=0;

		String transName;

		for(int i=0; i<nTargets; i++){
			nHits=0;
			nNewHits=0;

			br = new BufferedReader ( new InputStreamReader(new FileInputStream(   new File(annotationFile) ) ));
			tempS=br.readLine();

			while(tempS!=null){
				TEMPS=tempS.toUpperCase();
				int nPos=0;
				int nNeg=0;
				int nPosFound=0;
				int nNegFound=0;
				int minKeywordPos=99999;
				for(int j=0; j<nKeywords[i]; j++){
					String keyword=keywords[i][j];
					boolean isPos=true;
					if(keyword.startsWith("!")){
						isPos=false;
						nNeg++;
						keyword=keyword.substring(1);
					}else{
						nPos++;
					}
					if(TEMPS.indexOf(keyword)>=0){
						if(isPos){
							minKeywordPos=Math.min(minKeywordPos,TEMPS.indexOf(keyword));
							nPosFound++;
						}
					}else{
						if(!isPos){
							nNegFound++;
						}
					}
				}
				boolean hit=false;
				if(type[i]==0){	//AND
					if(nPosFound==nPos && nNegFound==nNeg){
						hit=true;
					}
				}else{			//OR
					if(nPosFound>0){
						hit=true;
					}
				}

				if(hit){
					//System.out.println(minKeywordPos);
					targetCounts[i]++;
					nHits++;
					nHitsTotal++;
					transName=tempS.split("\t")[0];
					candID=-1;
					for(int j=0; j<nExpected; j++){
						if(candidateNames[j]==null){
							break;
						}else if(candidateNames[j].equals(transName)){
							candID=j;
							break;
						}
					}

					if(candID<0){
						candidateKeywordIDs[nCandidates]=i;
						candidateKeywordPos[nCandidates]=minKeywordPos;
						targetCountsNew[i]++;
						nNewHits++;
						nNewHitsTotal++;
						candID=nCandidates;
						nCandidates++;
						candidateNames[candID]=transName;
						candidateAnnot[candID]=tempS;
					}
					candidateCounts[candID][i]++;
				}
				tempS=br.readLine();
			}
			System.out.println(nHits+"\t"+nNewHits+"\t"+targets[i]);
		}
		System.out.println(nHitsTotal+" total hits.");
		System.out.println(nNewHitsTotal+" unique hits.");

		BufferedWriter bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outStem+"_IDs.txt") ) ));	//out file
		for(int i=0; i<nCandidates; i++){
			bw.write(candidateNames[i]+"\n");
		}
		bw.flush();
		bw.close();

		bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outStem+"_TargetCounts.txt") ) ));	//out file
		bw.write("Count\tTarget\n");
		for(int i=0; i<nTargets; i++){
			bw.write(targetCounts[i]+"\t"+targetCountsNew[i]+"\t"+targets[i]+"\n");
		}
		bw.flush();
		bw.close();

		bw = new BufferedWriter ( new OutputStreamWriter(new FileOutputStream(   new File(outStem+"_Annotation.txt") ) ));	//out file
		for(int i=0; i<nCandidates; i++){
			for(int j=0; j<nTargets; j++){
				bw.write(candidateCounts[i][j]+"\t");				
			}
			bw.write(candidateKeywordIDs[i]+"\t"+candidateKeywordPos[i]+"\t"+candidateAnnot[i].substring(candidateKeywordPos[i],Math.min(candidateAnnot[i].length(),candidateKeywordPos[i]+50)));
			bw.write("\t"+candidateAnnot[i]+"\n");
		}
		bw.flush();
		bw.close();


      }catch(IOException ioe){System.out.println("<<!!ERROR main()!!>>"+ioe.getMessage());}
  }

}
