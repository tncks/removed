package modi;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class ProteinMatchInfo extends TreeMap<Integer, ArrayList<PeptideMatchInfo>> {

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

}
