package cn.hnu.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;
import java.util.Vector;

import cn.hnu.beans.Network;
import cn.hnu.beans.Node;

public class CFOCM_ATTACH {

	/**
	 * @param args
	 */
	 public static ArrayList<Node> proteins;
	public static  HashMap<String,Integer> pro_idx;
	 public static HashMap<String, Integer> goitems_idx; 
	public  static double[][] goitemssim;
	 
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		
		CFOCM_ATTACH alg = new CFOCM_ATTACH();
		
	}
	public static Vector<TreeSet<String>> Append_Core_Attachments(String corefile, String complexesfile)
	{
		proteins  = Network.proteins;
		pro_idx = Network.pro_idx;
		goitems_idx = Network.goitems_idx;
		goitemssim =Network.goitemssim;
		
		Vector<TreeSet<String>> complexes = GetCoreClqsFromFile(corefile);
		complexes = FindAttachmentsofCoreClqs(complexes);
		FileOperation.SaveComplexs2File(complexesfile, complexes);
		
		return complexes;
	
	}


	public static Vector<TreeSet<String>> FindAttachmentsofCoreClqs(Vector<TreeSet<String>> clqs)
	{
		for(int i=0;i<clqs.size();i++)

		{
			//System.out.println(i+" begins:");
			TreeSet<String> clq = clqs.get(i);
			int ix = FindOneAttachment(clq);
			while(ix!=-1)
			{
				clq.add(proteins.get(ix).getName());
				ix = FindOneAttachment(clq);
			}
//			unGOCCA
			if(!JudeOneClq(clq))
			{
				clqs.remove(i);
				i--;
			}
		}
		//clqs = NCAPM_CORE.DelDuplicatedClqs(clqs, 0.5);
		return clqs;
		
	}
	

	
	private static int FindOneAttachment(TreeSet<String> clq)
	{
		
		String[] pros = new String[clq.size()];
		clq.toArray(pros);
		
		int[] ipros = new int[pros.length];
		int[] addpros = new int[pros.length+1];
		for(int i=0;i<pros.length;i++)
		{
			ipros[i] = pro_idx.get(pros[i]);
			addpros[i] = ipros[i];
		}

		double orgGraphEntropy = CalculateFitness(ipros);
//		System.out.println(orgGraphEntropy);
		int mxAddix = -1;
		double mxGraphEntropy = orgGraphEntropy;
		double addGraphEntropy = 0;
		
		ArrayList<Integer> neighbors = FindNeighbors(clq);
		
//		double[] priorV = new double[neighbors.size()];
		
		for(int i=0;i<neighbors.size();i++)
		{
			int neighbor = neighbors.get(i);
			addpros[pros.length] = neighbor;
			addGraphEntropy = CalculateFitness(addpros);
		//	if(addGraphEntropy>mxGraphEntropy&&HasCommonGoFucs(clq,neighbor))
			if(addGraphEntropy>mxGraphEntropy)
			{
				mxGraphEntropy = addGraphEntropy;
				mxAddix = neighbor;
			}
		}
	//	System.out.println(" MX:"+mxGraphEntropy+" IDX:"+mxAddix);
//		System.out.println(mxAddix);
		return mxAddix;
	}
	
	public static boolean HasCommonGoFucs(TreeSet<String> clq, int addix)
	{
		ArrayList<String> mxCommonGoitem = findComplexCommonGoitems(clq);
		HashSet<String> golink = proteins.get(addix).getGolink();
//		System.out.println(addix+": ");
		if(golink.size()==0)
			return true;
		String[] sgolink = new String[golink.size()];
//		System.out.println(mxCommonGoitem.size()+"xxx"+golink.size());
		golink.toArray(sgolink);
		
		int[] igolink  = new int[sgolink.length];
		for(int i=0;i<igolink.length;i++)
			igolink[i] = goitems_idx.get(sgolink[i]);
		                 
		for(int i=0;i<mxCommonGoitem.size();i++)
		{
			int item = goitems_idx.get(mxCommonGoitem.get(i));
			if(golink.contains(mxCommonGoitem.get(i)))
				return true;
			else
			{
				for(int j=0;j<igolink.length;j++)
				{
					if(goitemssim[item][igolink[j]]>1.96)
						return true;
				}
			}
		}
		return false;
	}
	
	public static ArrayList<String> findComplexCommonGoitems(TreeSet<String> clq)
	{
		String[] pro = new String[clq.size()];
		clq.toArray(pro);
		
		ArrayList<Node> proteins = Network.proteins;
		HashMap<String,Integer> pro_idx = Network.pro_idx;
		
		HashSet<String> clqGOs = new HashSet<String>();
		
		int[] ipro = new int[pro.length];
		for(int i=0;i<pro.length;i++)
		{
			ipro[i] = pro_idx.get(pro[i]);
			HashSet<String> goitems=  proteins.get(ipro[i]).getGolink();
//			System.out.println(pro[i]+":"+goitems.size());
			clqGOs.addAll(goitems);
		}
		
		int[] count = new int[clqGOs.size()];
		String[] strclqGos = new String[clqGOs.size()];
		clqGOs.toArray(strclqGos);
//		System.out.println("ALL GO ITEMS:"+clqGOs.size());
		
		for(int i=0;i<pro.length;i++)
		{
			HashSet<String> goitems=  proteins.get(ipro[i]).getGolink();
			for(int j=0;j<clqGOs.size();j++)
			{
				if(goitems.contains(strclqGos[j]))
				{
					count[j]++;
				}
			}
		}
		
		int max = 0;
		for(int i=0;i<count.length;i++)
		{
//			System.out.println(count[i]);
			max = Math.max(count[i], max);
		}

		//System.out.println(max);
//    System.out.println("proteins num:"+pro.length);

		ArrayList<String> mxCommonGoitem = new ArrayList<String>();
		for(int i=0;i<count.length;i++)
		{
			if(count[i]==max)
			{
				mxCommonGoitem.add(strclqGos[i]);
			}
		}
		return mxCommonGoitem;
	}
	
	
	private static double CalculateFitness(int[] v)
	{
		double density = 0;
		double fitV = 0;
//		double[] prefix  = new double[v.length];
		for(int i=0;i<v.length;i++)
		{
			int[] neighbors2 = proteins.get(v[i]).getNeighbors();
			int degree = neighbors2.length;
			int indegree = intersect(v, neighbors2);
			density += indegree;
			double pi = indegree*1.0/degree;
			fitV += pi;
		}
		fitV = fitV*density/v.length/(v.length-1)/v.length;
		return fitV;
	}
	
	

	private static int intersect(int[] a, int[] b)
	{
		int common = 0;
		for(int i=0;i<a.length;i++)
			for(int j=0;j<b.length;j++)
			{
				if(a[i]==b[j])
					common++;
			}
		return common;
	}
	
	private static ArrayList<Integer> FindNeighbors(TreeSet<String> clq)
	{
		ArrayList<Integer> neighbors = new ArrayList<Integer>();
		Iterator<String> it = clq.iterator();
		while(it.hasNext())
		{
			String pro = it.next();
			int ix = pro_idx.get(pro);
			int[] neighs = proteins.get(ix).getNeighbors();
			
			for(int i=0;i<neighs.length;i++)
				if(!neighbors.contains(neighs[i])&&!clq.contains(proteins.get(neighs[i]).getName()))
					neighbors.add(neighs[i]);
		}
		return neighbors;
	}
	
	
	private static Vector<TreeSet<String>> GetCoreClqsFromFile(String filename)
	{
		Vector<TreeSet<String>> clqs = new Vector<TreeSet<String>>();
		try {
			FileReader f = new FileReader(filename);
			BufferedReader in = new BufferedReader(f);
			String s = in.readLine();
			while (s != null) {
				String[] temp = s.split(" ");
				TreeSet<String> clq = new TreeSet<String>();
				for (int i = 0; i < temp.length; i++) {
					clq.add(temp[i]);
				}
				clqs.add(clq);
				s = in.readLine();
			}
		} catch (Exception e) {
		}
		return clqs;
	}
	
	private static boolean JudeOneClq(TreeSet<String> clq)
	{
		String[] pro = new String[clq.size()];
		clq.toArray(pro);
		
		ArrayList<Node> proteins = Network.proteins;
		HashMap<String,Integer> pro_idx = Network.pro_idx;
		
		HashSet<String> clqGOs = new HashSet<String>();
		
		int[] ipro = new int[pro.length];
		for(int i=0;i<pro.length;i++)
		{
			ipro[i] = pro_idx.get(pro[i]);
			HashSet<String> goitems=  proteins.get(ipro[i]).getGolink();
//			System.out.println(pro[i]+":"+goitems.size());
			clqGOs.addAll(goitems);
		}
		
		int[] count = new int[clqGOs.size()];
		String[] strclqGos = new String[clqGOs.size()];
		clqGOs.toArray(strclqGos);
//		System.out.println("ALL GO ITEMS:"+clqGOs.size());
		
		for(int i=0;i<pro.length;i++)
		{
			HashSet<String> goitems=  proteins.get(ipro[i]).getGolink();
			for(int j=0;j<clqGOs.size();j++)
			{
				if(goitems.contains(strclqGos[j]))
				{
					count[j]++;
				}
			}
		}
		
		int max = 0;
		for(int i=0;i<count.length;i++)
		{
//			System.out.println(count[i]);
			max = Math.max(count[i], max);
		}
//System.out.println("max"+max);
    	if(max <= clq.size()*1.0/2)
    		return false;
    	else
    		return true;
	}
	
}
