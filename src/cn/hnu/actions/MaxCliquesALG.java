package cn.hnu.actions;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.TreeSet;
import java.util.Vector;

import cn.hnu.beans.Network;
import cn.hnu.beans.Node;


public class MaxCliquesALG {

	/**
	 * @param args
	 */
	public static List<Vector<TreeSet<String>>> GetMaxCliques(boolean isDIP, boolean isFirst){
		
		String prefix ; 
		String mxprefix; 
		
		int mxclqsize = 0;
		if(isDIP)
		{
			prefix = "cliques/dip_clqs";
			mxprefix = "cliques/dip_mxclqs";
			mxclqsize = 10;
		}
		else
		{
			mxclqsize = 12;
			prefix = "cliques/gavin_clqs";
			mxprefix = "cliques/gavin_mxclqs";
		}
		
		List<Vector<TreeSet<String>>> clqs = new ArrayList<Vector<TreeSet<String>>>();
		
		if(!isFirst){
			for (int i = 3; i <= mxclqsize; i++) {
				String path = mxprefix + String.valueOf(i);
				clqs.add(FileOperation.ReadComplexsFromFile(path));
			}
			return clqs;
		}
		
//		FindCliquesOfSize(prefix,
//				Network.graph,
//				Network.proteins, Network.pro_idx, mxclqsize);
		
		for (int i = 3; i <= mxclqsize; i++) {
			String path = prefix + String.valueOf(i);
			clqs.add(FileOperation.ReadComplexsFromFile(path));
		}
		
		// /删除非最大团
		Vector<TreeSet<String>> higherclqs = clqs.get(clqs.size() - 1);
		String filename = mxprefix + String.valueOf(clqs.size()-1+3);
		
		FileOperation.SaveComplexs2File(filename, higherclqs);
		
		for (int i = clqs.size() - 2; i >= 0; i--) {
			Vector<TreeSet<String>> curclq = clqs.get(i);
			curclq = DelDuplicatedClqs(curclq, higherclqs);
			clqs.set(i, curclq);
			FileOperation.SaveComplexs2File(mxprefix+String.valueOf(i+3), curclq);
			higherclqs = curclq;
//			System.out.println("Yes");
		}
		
		return clqs;
		
	}
	

	private static Vector<TreeSet<String>> FindCliquesOfSize(int start, int end,
			String pathprefix,
			int[][] graph,
			ArrayList<Node> proteins, HashMap<String, Integer> pro_idx, int size) {

		if (size == 2) {
			String filename = pathprefix + String.valueOf(size);
			Vector<TreeSet<String>> clqs2 = FileOperation.ReadComplexsFromFile(filename);
			return clqs2;
		}

		Vector<TreeSet<String>> clqs = new Vector<TreeSet<String>>();
		HashSet<TreeSet<String>> hset = new HashSet<TreeSet<String>>();
		if (size >= 3) {
			String filename = pathprefix + String.valueOf(size-1);
			Vector<TreeSet<String>> mclqs = FileOperation.ReadComplexsFromFile(filename);
			System.out.println(mclqs.size());
			
			for (int i = start; i < end; i++) {
				TreeSet<String> ni = mclqs.elementAt(i);

				System.out.println(i);
				for (int j = i + 1; j < mclqs.size(); j++) {
					// /ni-nj

					TreeSet<String> set = new TreeSet<String>();
					TreeSet<String> nj = mclqs.elementAt(j);
					set.clear();
					set.addAll(ni);
					set.removeAll(nj);

					if (set.size() != 1)
						continue;

					String ti = set.first();
					// System.out.print("here1: "+ti+" "+set.size());
					// nj-ni
					set = new TreeSet<String>();
					set.addAll(nj);
					set.removeAll(ni);

					if (set.size() != 1)
						continue;
					String tj = set.first();

					// System.out.println("HERE2: "+ti+"-"+tj);

					// System.out.println("HERE2: "+ti+"-"+tj);

					if(pro_idx.containsKey(ti)&&pro_idx.containsKey(tj))
					if (graph[pro_idx.get(ti)][pro_idx.get(tj)] == 1) { // ni U
						// nj
						//System.out.println("HERE3");
						set = new TreeSet<String>();
						set.addAll(ni);
						set.addAll(nj);
						hset.add(set);
					//	if (!clqs.contains(set)) {
						//	clqs.add(set);
							// System.out.println(set);
						//}

					}
				}
			}
		}
		clqs.addAll(hset);
		String filename = pathprefix + String.valueOf(size)+"_"+String.valueOf(start)+"_"+String.valueOf(end);
		FileOperation.SaveComplexs2File(filename, clqs);
		return clqs;
	}
	
	
	
	private static Vector<TreeSet<String>> FindCliquesOfSize(
			String pathprefix,
			int[][] graph,
			ArrayList<Node> proteins, HashMap<String, Integer> pro_idx, int size) {

		if (size == 2) {
			Vector<TreeSet<String>> clqs2 = new Vector<TreeSet<String>>();
			for (int i = 0; i < graph.length; i++) {
				for (int j = i + 1; j < graph[0].length; j++) {
					TreeSet<String> e = new TreeSet<String>();
					if (graph[i][j] == 1) {
						e.add(proteins.get(i).getName());
						e.add(proteins.get(j).getName());
						clqs2.addElement(e);
					}

				}
			}
			String filename = pathprefix + String.valueOf(size);
			FileOperation.SaveComplexs2File(filename, clqs2);
			return clqs2;
		}

		Vector<TreeSet<String>> clqs = new Vector<TreeSet<String>>();
		HashSet<TreeSet<String>> hset = new HashSet<TreeSet<String>>();
		if (size >= 3) {
			Vector<TreeSet<String>> mclqs = FindCliquesOfSize(pathprefix, graph, proteins,
					pro_idx, size - 1);
			//System.out.println("22222:" + mclqs.size());
			for (int i = 0; i < mclqs.size(); i++) {
				TreeSet<String> ni = mclqs.elementAt(i);

				System.out.println(i);
				for (int j = i + 1; j < mclqs.size(); j++) {
					// /ni-nj

					TreeSet<String> set = new TreeSet<String>();
					TreeSet<String> nj = mclqs.elementAt(j);
					set.clear();
					set.addAll(ni);
					set.removeAll(nj);

					if (set.size() != 1)
						continue;

					String ti = set.first();
					// System.out.print("here1: "+ti+" "+set.size());
					// nj-ni
					set = new TreeSet<String>();
					set.addAll(nj);
					set.removeAll(ni);

					if (set.size() != 1)
						continue;
					String tj = set.first();

					// System.out.println("HERE2: "+ti+"-"+tj);

					if(pro_idx.containsKey(ti)&&pro_idx.containsKey(tj))
					if (graph[pro_idx.get(ti)][pro_idx.get(tj)] == 1) { // ni U
						// nj
						//System.out.println("HERE3");
						set = new TreeSet<String>();
						set.addAll(ni);
						set.addAll(nj);
						hset.add(set);
					//	if (!clqs.contains(set)) {
						//	clqs.add(set);
							// System.out.println(set);
						//}

					}
				}
			}
		}
		clqs.addAll(hset);
		String filename = pathprefix + String.valueOf(size);
		FileOperation.SaveComplexs2File(filename, clqs);
		return clqs;
	}

	
	public static Vector<TreeSet<String>> DelDuplicatedClqs(
			Vector<TreeSet<String>> curclqs, Vector<TreeSet<String>> higherclqs) {
		for (int i = 0; i < curclqs.size(); i++) {
			boolean flag = false;
			for (int j = 0; j < higherclqs.size(); j++) {
				if (higherclqs.get(j).containsAll(curclqs.get(i))) {
					flag = true;
					break;
				}

			}
			if (flag) {
				curclqs.remove(i);
				i--;
			}
		}
		return curclqs;
	}
	
	
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub

		PPNHandler.init_ppn("ppn/Gavin.txt");
		List<Vector<TreeSet<String>>> clqs = GetMaxCliques(false, true);
		// 对网络中的蛋白质进行GO注释, 产生蛋白质的GOlink,GO item到 id编号的映射:goitems_idx
		//GoAnnotation.Pros_GoAnnotaion("go/go_note.txt");
//		FindCliquesOfSize(200000,250000,
//				"cliques/Gavin_clqs",
//				Network.graph,
//				Network.proteins, Network.pro_idx, 5);
//		Vector<TreeSet<String>> clqs = FileOperation.ReadComplexsFromFile("cliques/biogrid_clqs4");
//		System.out.println(clqs.size());
//		FileOperation.SaveComplexs2File("cliques/biogrid_clq42", clqs);
	}

}
