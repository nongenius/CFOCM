package cn.hnu.actions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import cn.hnu.beans.Network;
import cn.hnu.beans.Node;




/**
 * @author libo
 *
 *	backup 2013.10.2
 */
public class PPNHandler {

//	public static HashMap<String,Integer> pro_idx; //蛋白质到编号id
//	public static HashMap<Integer,String> idx_pro; //id到蛋白质
	
//	public static int[][] graph;
	
	public static void main(String args[]) throws Exception
	{
		PPNHandler.init_ppn("ppn/dip.txt");
	}
	

	
	public static void init_ppn(String ppnfile) throws Exception
	{
		///the number of proteins in DIP or Gavin
		Network.ProteinNUM = 4930;
//		Network.ProteinNUM = 1430;
		ArrayList<Node> proteins = new ArrayList<Node>();
		int[][] graph = new int[Network.ProteinNUM][Network.ProteinNUM]; ///dip中包括4930个蛋白质
		HashMap<String,Integer> pro_idx = new HashMap<String,Integer>();
		
		BufferedReader br = new BufferedReader(new FileReader(ppnfile));
		String sin;
		int idx = 0;
		int PPINUM =0;
		while((sin = br.readLine())!=null)
		{
			PPINUM++;
			String[] temp = sin.split(" |\t");
			if(!pro_idx.containsKey(temp[0]))
			{
				Node n = new Node(temp[0],idx);
				proteins.add(n);
				pro_idx.put(temp[0], idx);
				idx++;
			}
			if(!pro_idx.containsKey(temp[1]))
			{
				Node n = new Node(temp[1],idx);
				proteins.add(n);
				pro_idx.put(temp[1], idx);
				idx++;
			}
			graph[pro_idx.get(temp[0])][pro_idx.get(temp[1])] = 1;
			graph[pro_idx.get(temp[1])][pro_idx.get(temp[0])] = 1;	
		}
		br.close();
		
		Network.graph = graph;
		Network.proteins = proteins;
		Network.pro_idx = pro_idx;
		Network.PPIsNUM = PPINUM;
		System.out.println(Network.PPIsNUM);
		
		set_proteins_neighbors();
	}
	
	
	public static void set_proteins_neighbors()
	{
		ArrayList<Node> proteins = Network.proteins; 
		int[][] graph = Network.graph;
		for(int i=0;i<graph.length;i++)
			for(int j=i+1;j<graph[0].length;j++)
			{
				if(graph[i][j]==1)
				{
					int[] neighborsi = proteins.get(i).getNeighbors();
					if(neighborsi==null)
					{
						neighborsi = new int[1];
						neighborsi[0] = j;
					}
					else
					{
						int[] temp = new int[neighborsi.length+1];
						for(int k=0;k<neighborsi.length;k++)
							temp[k] = neighborsi[k];
						temp[neighborsi.length] = j;
						neighborsi = temp;
					}
					proteins.get(i).setNeighbors(neighborsi);
					proteins.get(i).setDegree(proteins.get(i).getDegree()+1);
					
					int[] neighborsj = proteins.get(j).getNeighbors();
					if(neighborsj==null)
					{
						neighborsj = new int[1];
						neighborsj[0] = i;
					}
					else
					{
						int[] temp = new int[neighborsj.length+1];
						for(int k=0;k<neighborsj.length;k++)
							temp[k] = neighborsj[k];
						temp[neighborsj.length] = i;
						neighborsj = temp;
					}
					proteins.get(j).setNeighbors(neighborsj);
					proteins.get(j).setDegree(proteins.get(j).getDegree()+1);
				}
			}
	}
	
	
	public static void init_ppn2(String ppnfile) throws Exception
	{
//		pro_idx = new HashMap<String,Integer>();
//		idx_pro = new HashMap<Integer,String>();
		
		HashMap<String,Integer> pro_idx = new HashMap<String,Integer>();
		HashMap<Integer,Node> idx_node = new HashMap<Integer,Node>();
		
		int[][] graph = new int[4930][4930]; ///dip中包括4930个蛋白质
		
		BufferedReader br = new BufferedReader(new FileReader(ppnfile));
		String sin;
		int idx = 0;
		while((sin = br.readLine())!=null)
		{
			String[] temp = sin.split(" |\t");
			if(!pro_idx.containsKey(temp[0]))
			{
				Node n = new Node(temp[0],idx);
				idx_node.put(idx, n);
				pro_idx.put(temp[0], idx);
//				idx_pro.put(idx, temp[0]);
				idx++;
			}
			if(!pro_idx.containsKey(temp[1]))
			{
				Node n = new Node(temp[1],idx);
				idx_node.put(idx, n);
				pro_idx.put(temp[1], idx);
//				idx_pro.put(idx, temp[1]);
				idx++;
			}
			graph[pro_idx.get(temp[0])][pro_idx.get(temp[1])] = 1;
			graph[pro_idx.get(temp[1])][pro_idx.get(temp[0])] = 1;	
		}
		br.close();
		
//		Network.idx_node = ;
		Network.graph = graph;
		Network.pro_idx = pro_idx;
	}
}

