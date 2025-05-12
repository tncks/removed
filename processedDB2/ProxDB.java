package processedDB;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import modi.Constants;
import modi.Mutables;

public class ProxDB extends ArrayList<Prox> {
	
	int sizeOfEntries;
	int sizeOfResidues;	
	String name, orgName;
	
	public int getSizeOfResidues() { return sizeOfResidues; }
	public String getName() { return orgName; }
	
	public void setSizeOfEntries(int x) { sizeOfEntries = x; }
	public void setSizeOfResidues(int x) { sizeOfResidues = x; }
	
	public void readFasta(String fileName) {
		if( fileName == null ) {
			System.out.println( "The protein DB was not specified." );
			return;
		}
		orgName = fileName;
		int residues= 0;
		try 
		{		
			if( fileName.lastIndexOf('.') > -1 )
				name= fileName.substring(0, fileName.lastIndexOf('.'));
			else name = fileName;
			
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			String s= "<";
			StringBuffer buffer;	
			
			while( s.startsWith(">") == false ) {
				s = in.readLine();
				if( s == null ) break;
			}
		
			while( s != null ) 
			{
				Prox protein= new Prox();				
				protein.setHeader( s.substring(1) );			
				
				buffer = new StringBuffer();
				while( (s = in.readLine()) != null ){
					if( s.startsWith(">") ) break;

					for( int aa=0; aa<s.length(); aa++ ){
						if( Character.isLetter( s.charAt(aa) ) )
							buffer.append( Character.toUpperCase(s.charAt(aa)) );
					}
				}
				if( buffer.length() < 3 ) continue; // check sequence
				
				protein.setSequence( buffer.toString() );
				residues += buffer.length();
				this.add( protein );					
			}
			in.close();						
		} 
		catch (FileNotFoundException e) 
		{
			System.out.println( "Cannot find the protein file, "+fileName );
			e.printStackTrace();
			System.exit(1);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		sizeOfEntries = this.size();
		sizeOfResidues = residues;	
		System.out.println( sizeOfEntries+" proteins / "  + sizeOfResidues+" residues" );
		
		if( Constants.targetDecoy == 1 ){
			System.out.print("Decoy search checked.....  ");

			for(int i=0; i<sizeOfEntries; i++){
				this.add( new Prox(Constants.DECOY_LABEL+this.get(i).getAccession(), this.get(i).getReverseSequence()) );
			}
			System.out.println( "generated " + sizeOfEntries+" decoy proteins");
			sizeOfEntries *= 2;
			sizeOfResidues *= 2;	
		}
	}
	
	public int getProteinIndex(String accession) {
		int index = -1;
		for(int i=0; i<this.size(); i++){
			if( accession.compareTo(get(i).getAccession()) == 0 ){
				index = i;
				break;
			}
		}
		return index;
	}
}























