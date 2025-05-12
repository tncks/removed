package modi;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class ProteinMatchInfo extends TreeMap<Integer, ArrayList<PeptideMatchInfo>> {
	public void addPeptideInfo(int proteinID, PeptideMatchInfo pepMatchInfo)
	{
		ArrayList<PeptideMatchInfo> pepMatchList = this.get(Integer.valueOf(proteinID));
		if(pepMatchList == null)
		{
			pepMatchList = new ArrayList<PeptideMatchInfo>();
			pepMatchList.add(pepMatchInfo);
			put(Integer.valueOf(proteinID), pepMatchList);
		}
		else
		{
			pepMatchList.add(pepMatchInfo);
		}
	}
	
	public String toString()
	{
		StringBuffer output = new StringBuffer();
		output.append("Protein Match Summary\n");
		Iterator<Map.Entry<Integer, ArrayList<PeptideMatchInfo>>> it = (this.entrySet()).iterator();
		while(it.hasNext())
		{
			Map.Entry<Integer, ArrayList<PeptideMatchInfo>> entry = it.next();
			output.append("Protein : " + entry.getKey() + "\n");
			output.append(entry.getValue() + "\n");
		}
		return output.toString();
	}
	
	public ProteinContainer getMatchedProtein(ProteinContainer proteins)
	{
		ProteinContainer matchedProtein = new ProteinContainer();
		Iterator<Integer> it = this.keySet().iterator();
		while(it.hasNext()) {
			int index = it.next().intValue();
			matchedProtein.add(proteins.get(index));
		}
		return matchedProtein;
	}
}
