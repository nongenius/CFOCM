package cn.hnu.actions;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import cn.hnu.beans.Network;
import cn.hnu.beans.Node;



public class GoAnnotation {

	/**
	 * @param libo
	 */
	



	public static double[][] Compute_GoItems_Similar()
	{
		int EdgeNUM = Network.PPIsNUM;
		System.out.println(EdgeNUM);
		//		double obs,exp,r;
		double[][] obs = new double[Network.GoitemsNUM][Network.GoitemsNUM];
		double[][] exp = new double[Network.GoitemsNUM][Network.GoitemsNUM];

		double[][] goitemssim = new double[Network.GoitemsNUM][Network.GoitemsNUM];

		for(int i=0;i<goitemssim.length;i++)
			for(int j=0;j<goitemssim[0].length;j++)
			{
				goitemssim[i][j]=-1;
			}

		HashMap<String, HashSet<Integer>> goitem_pros = Network.goitem_pros;
		HashMap<String, Integer> goitems_idx = Network.goitems_idx ;
		ArrayList<Node> proteins = Network.proteins;
		int graph[][]= Network.graph;

		Set<String> set = goitem_pros.keySet();
		String[] itemset = new String[set.size()];
		set.toArray(itemset);
		
		int item1,item2;
		int n1,n2;
		for(int i=0;i<itemset.length;i++)
//		for(int i=0;i<1;i++)
		{
			item1 = goitems_idx.get(itemset[i]);
//			System.out.println("item1:"+item1);
//			System.out.println(i+" "+item1);
			HashSet<Integer> ps1 = goitem_pros.get(itemset[i]);
//			System.out.println("size1:"+ ps1.size());
//			n1 = ps1.size();
//			Integer[] p1 = (Integer[]) ps1.toArray();
			Integer[] p1 = new Integer[ps1.size()];
			ps1.toArray(p1);
			for(int j=i;j<itemset.length;j++)
//			for(int j=i+1;j<2;j++)
			{
				item2 = goitems_idx.get(itemset[j]);
				HashSet<Integer> ps2 = goitem_pros.get(itemset[j]);
				
				
//				System.out.println("size2:"+ ps2.size());
//				n2 = ps2.size();
				Integer[] p2 = new Integer[ps2.size()];
				ps2.toArray(p2);
				n1=n2=0;
				for(int k=0;k<p1.length;k++)
				{
					n1+=proteins.get(p1[k]).getDegree();
				}
				for(int l=0;l<p2.length;l++)
				{
					n2+=proteins.get(p2[l]).getDegree();
				}
				for(int k=0;k<p1.length;k++)
					for(int l=0;l<p2.length;l++)
					{
							if(graph[p1[k]][p2[l]]==1)
							{
								obs[item1][item2]++;
								obs[item2][item1]++;
							}
					}
				
//				System.out.println("degree sum"+ n1+": "+n2);
//				System.out.println("item2:"+item2);
				exp[item1][item2] = n1*n2*1.0/EdgeNUM;
//				if(exp[item1][item2]>1)
//					System.out.println("XXX"+item1+" "+item2+" "+exp[item1][item2]);
				exp[item2][item1] = exp[item1][item2];
			}
		}

		
		double[] sobs = new double[Network.GoitemsNUM];
		for(int i=0;i<Network.GoitemsNUM;i++)
//		for(int i=0;i<2;i++)
		{
			sobs[i] = 0;
			
			for(int k=0;k<Network.GoitemsNUM;k++)
			{
				sobs[i] += obs[i][k];
			}
			
		}

		int count =0;
		int count2= 0;
		for(int i=0;i<Network.GoitemsNUM;i++)
			for(int j=i;j<Network.GoitemsNUM;j++)
			{
				if(obs[i][j]>exp[i][j])
					goitemssim[i][j] = (obs[i][j]-exp[i][j])/Math.sqrt(Math.abs(exp[i][j]*(1-sobs[i]/EdgeNUM)*(1-sobs[j]/EdgeNUM)));
//				
				if(goitemssim[i][j]!=-1)
				{
					count++;
					System.out.println(count+" "+goitemssim[i][j]);
				}
				else
				{
					count2++;
				}
				
				if(Double.isNaN(goitemssim[i][j]))
				{
					goitemssim[i][j]=-1;
						System.out.println(goitemssim[i][j]+"-"+i+"-"+j+"-"+obs[i][j]+"-"+exp[i][j]+"-"+sobs[i]+"-"+sobs[j]);
				}
						goitemssim[j][i] = goitemssim[i][j];
			}
		FileOperation.SaveMatrix2File(goitemssim,"dip_gosim.txt");
		return goitemssim;
	}


	private static ArrayList<Integer> intersect(int[] a, int[] b)
	{
		ArrayList<Integer> c = new ArrayList<Integer>();
		for(int i=0;i<a.length;i++)
			for(int j=0;j<b.length;j++)
			{
				if(a[i]==b[j])
					c.add(a[i]);
			}
		return c;
	}
	
	
	public static void Pros_GoAnnotaion(String filePath) throws Exception  
	{
		int GoitemsNUM = 0;

		HashMap<String,Integer> pro_idx = Network.pro_idx;
		ArrayList<Node> proteins = Network.proteins;
		HashMap<String, Integer> goitems_idx =  new HashMap<String, Integer>();

		HashMap<String, HashSet<Integer>> goitem_pros = new HashMap<String, HashSet<Integer>>();

		BufferedReader bw = new BufferedReader(new FileReader(filePath));
		String temp;
		int idx;

		while((temp = bw.readLine()) != null)
		{
			String[] edges = temp.split(" ");
			if(pro_idx.containsKey(edges[0]))
			{
				idx = pro_idx.get(edges[0]);
				Node n = proteins.get(idx);
				HashSet<String> goa =  n.getGolink();
				goa.add(edges[1]);
				n.setGolink(goa);
				if(goitem_pros.containsKey(edges[1]))
				{
					HashSet<Integer> pros = goitem_pros.get(edges[1]);
					pros.add(idx);
				}
				else
				{
					HashSet<Integer> pros = new HashSet<Integer>(); 
					pros.add(idx);
					goitem_pros.put(edges[1], pros);
				}
				if(!goitems_idx.containsKey(edges[1]))
				{
					goitems_idx.put(edges[1], GoitemsNUM);
					GoitemsNUM++;
				}
			}
		}
		bw.close();

		Network.GoitemsNUM = GoitemsNUM;
		Network.goitems_idx = goitems_idx;
		Network.goitem_pros = goitem_pros;

	}
}

