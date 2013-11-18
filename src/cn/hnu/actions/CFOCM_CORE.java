package cn.hnu.actions;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;
import java.util.Vector;

import cn.hnu.beans.Network;
import cn.hnu.beans.Node;

public class CFOCM_CORE {

	/**
	 * @param args
	 * @throws Exception
	 */

	static Vector<TreeSet<String>> cclqs = new Vector<TreeSet<String>>();
	private static Vector<TreeSet<String>> higherclqs;

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		
		PPNHandler.init_ppn("ppn/dip.txt");
		// 对网络中的蛋白质进行GO注释, 产生蛋白质的GOlink,GO item到 id编号的映射:goitems_idx
		GoAnnotation.Pros_GoAnnotaion("go/go_note.txt");
		Network.goitemssim = FileOperation.ReadMatrixFromFile("gosim.txt");
		//Network.pro_gosim = FileOperation.ReadMatrixFromFile("pro_gosim.txt");
		CFOCM_CORE.SearchTheBestMergeableThres();
	}

	public static void test() {
		Vector<TreeSet<String>> clqs = FileOperation
				.ReadComplexsFromFile("test.txt");
		System.out.println(clqs.size());
		higherclqs = new Vector<TreeSet<String>>();
		// clqs = CombineCLQs3(clqs,0);
		System.out.println(clqs.size());
		System.out.println(clqs.get(0).size());
	}
	
	

	public static Vector<TreeSet<String>> Generate(int[][] graph,
			ArrayList<Node> proteins, HashMap<String, Integer> pro_idx, int size) {
		String path = "cliques/clqs";
		Vector<TreeSet<String>> clqs = null;
		if (new File(path).exists())// 如果已经生成clqs，读取即可
		{
			clqs = FileOperation.ReadComplexsFromFile(path);
		} else // 没有
		{
			if (!new File("cliques/clqs10").exists())// 不同size的团不存在，生成之
			{
				//CFOCM_CORE.FindCliquesOfSize(graph, proteins, pro_idx, 10);
			}
			// clqs = Detect_Complex_Cores(); //生成最终的种子团
		}
		return clqs;
	}

	public static void SearchTheBestMergeableThres() {
		
		String prefix = "cores/dip_cores_";
		String suffix = ".txt";

		// 检测cores
		double threshold = 0.2;
		do{
			List<Vector<TreeSet<String>>> clqs = MaxCliquesALG.GetMaxCliques(true, false);
			Vector<TreeSet<String>> cclqs = Detect_Complex_Cores(clqs, threshold);
			// 保存cores
			String savepath = prefix + String.valueOf(threshold) + suffix;
			FileOperation.SaveComplexs2File(savepath, cclqs);
			System.out.println("Cores saved to Path:" + savepath);
			threshold += 0.01;
			String sthreshold = String.format("%.2f",threshold);
			threshold = Double.parseDouble(sthreshold);
		}while(threshold<=0.5);
	}

	public static Vector<TreeSet<String>> Detect_Complex_Cores(
			List<Vector<TreeSet<String>>> clqs, double threshold) {
		
		
		System.out
				.println("Protein Complex Cores Detection Start, With Mergable threshold:"
						+ threshold);

		Vector<TreeSet<String>> cclqs = new Vector<TreeSet<String>>();
		higherclqs = new Vector<TreeSet<String>>();
		for (int i = clqs.size() - 1; i >= 0; i--) {
			Vector<TreeSet<String>> curclq = clqs.get(i);
			curclq = DelDuplicatedClqs(curclq, higherclqs);
			// System.out.println((i+3)+"Before:"+curclq.size());
			curclq = Merge_Similar_Cores(curclq, threshold);
			// System.out.println((i+3)+"After:"+curclq.size());
			clqs.set(i, curclq);
			cclqs.addAll(curclq);
			//higherclqs = curclq;
		}
		System.out.println(cclqs.size());
		

//		cclqs = Merge_Similar_Cores(cclqs, threshold);

		cclqs = DelDuplicatedClqs(cclqs, 0.5);
		
		
		System.out.println(cclqs.size());
//		FileOperation.SaveComplexs2File("cores/dip_cores_no_"+String.valueOf(threshold)+ ".txt",
//				cclqs);
		
		System.out
				.println("Protein Complex Cores Detection End, Returned Cores:"
						+ cclqs.size());
		return cclqs;
		
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

	public static Vector<TreeSet<String>> DelDuplicatedClqs(
			Vector<TreeSet<String>> clqs, double threshold) {
		for(int i=0;i<clqs.size();i++)
		{
			TreeSet<String> clq = clqs.get(i);
			if(clq.size()<3)
			{
				clqs.remove(i);
				i--;
			}
		}
		for(int i=0; i<clqs.size(); i++)
		{
			TreeSet<String> clq = clqs.get(i);
			for(int j=0; j<clqs.size(); j++)
			{
				if(i!=j&&clq.containsAll(clqs.get(j)))
				{
					clqs.remove(j);
					j--;
				}
			}
		}
		return clqs;
	}

	public static Vector<TreeSet<String>> Merge_Similar_Cores(
			Vector<TreeSet<String>> clqs, double threshold) {
		HashMap<String, Integer> pro_idx = Network.pro_idx;
		ArrayList<Node> proteins = Network.proteins;
		HashMap<String, Integer> goitems_idx = Network.goitems_idx;
		double[][] goitemssim = Network.goitemssim;
		Vector<TreeSet<String>> clqs_fiter;

		int prev,cur;
		//do{
			prev = clqs.size();
			Vector<TreeSet<Integer>> familys = new Vector<TreeSet<Integer>>();
			for (int i = 0; i < clqs.size(); i++) {
				TreeSet<Integer> family = new TreeSet<Integer>();
				family.add(i);
				for (int j = 0; j < clqs.size() && i >= 0; j++) {
					if(i==j)
						continue;
					TreeSet<String> set = new TreeSet<String>();
					set.addAll(clqs.get(i));
					set.retainAll(clqs.get(j));
					TreeSet<String> intersect = set;
					int common = intersect.size();
					int sizei = clqs.get(i).size();
					int sizej = clqs.get(j).size();
					int maxsize = Math.min(sizei, sizej);
					///sizei == sizej
					double similar = common*1.0*common/sizei/sizej;
					if(similar>=threshold)
						family.add(j);
				}
				//			System.out.println(i);
				//familys.insertElementAt(family, i);
				familys.add(family);
			}

			clqs_fiter= new Vector<TreeSet<String>>();

			for(int i=0;i<familys.size();i++)
			{
				boolean flag = false;
				for(int j=i+1;j<familys.size();j++)
					if(HasIntersection(familys.get(i),familys.get(j)))
					{
						familys.get(i).addAll(familys.get(j));
						familys.remove(j);
						j--;
						flag = true;
					}
				if(flag)
					i--;
			}

			for(int i=0;i<familys.size();i++)
			{
				TreeSet<Integer> family = familys.get(i);
				
				TreeSet<String> pros = new TreeSet<String>();
				Iterator<Integer> it = family.iterator();
				while(it.hasNext()){
					pros.addAll(clqs.get(it.next()));
				}

				TreeSet<String> pros_common = new TreeSet<String>();
				Iterator<String> it2 = pros.iterator();
				while(it2.hasNext())
				{
					String curpro = it2.next();
					int proidx = pro_idx.get(curpro);
					if (HasCommonGoFucs(pros, proidx, proteins, pro_idx,
							goitems_idx, goitemssim))
						pros_common.add(curpro);
				}
				
				int del = DeleteONE(pros_common, proteins, pro_idx);
				while(del!=-1&&pros_common.size()>2)
				{
					pros_common.remove(proteins.get(del).getName());
					del = DeleteONE(pros_common, proteins, pro_idx);
				}
				
				if(pros_common.size()>=3){
					
					clqs_fiter.add(pros_common);
					higherclqs.add(pros_common);
				}
				
			}
			cur = clqs_fiter.size();
			clqs = clqs_fiter;
		//}while(prev!=cur);
		
		return clqs_fiter;
	}

	public static boolean HasCommonGoFucs(TreeSet<String> clq, int addix,
			ArrayList<Node> proteins, HashMap<String, Integer> pro_idx,
			HashMap<String, Integer> goitems_idx, double[][] goitemssim) {

		// 找出COMMON GO ITEMS
		ArrayList<String> mxCommonGoitem = FindComplexCommonGoitems(clq,
				proteins, pro_idx);
		// 当前节点对应的GOLINK
		HashSet<String> golink = proteins.get(addix).getGolink();
		// System.out.println(addix+": ");
		// 如果没有注释
		if (golink.size() == 0)
			return false;

		// 转化成字符串数组
		String[] sgolink = new String[golink.size()];
		// System.out.println(mxCommonGoitem.size()+"xxx"+golink.size());
		golink.toArray(sgolink);

		// 映射到GO link id
		int[] igolink = new int[sgolink.length];
		for (int i = 0; i < igolink.length; i++)
			igolink[i] = goitems_idx.get(sgolink[i]);

		// 判断是否有共同的COMMON GO ITEM
		for (int i = 0; i < mxCommonGoitem.size(); i++) {
			int item = goitems_idx.get(mxCommonGoitem.get(i));
			if (golink.contains(mxCommonGoitem.get(i)))
				return true;
			else {
				for (int j = 0; j < igolink.length; j++) {
					if (goitemssim[item][igolink[j]] > 1.96)
						return true;
				}
			}
		}
		return false;
	}

	public static ArrayList<String> FindComplexCommonGoitems(
			TreeSet<String> clq, ArrayList<Node> proteins,
			HashMap<String, Integer> pro_idx) {
		String[] pro = new String[clq.size()];
		clq.toArray(pro);

		HashSet<String> clqGOs = new HashSet<String>();

		int[] ipro = new int[pro.length];
		for (int i = 0; i < pro.length; i++) {
			//System.out.println(pro[i]);
			ipro[i] = pro_idx.get(pro[i]);
			HashSet<String> goitems = proteins.get(ipro[i]).getGolink();
			// System.out.println(pro[i]+":"+goitems.size());
			clqGOs.addAll(goitems);
		}

		int[] count = new int[clqGOs.size()];
		String[] strclqGos = new String[clqGOs.size()];
		clqGOs.toArray(strclqGos);
		// System.out.println("ALL GO ITEMS:"+clqGOs.size());

		for (int i = 0; i < pro.length; i++) {
			HashSet<String> goitems = proteins.get(ipro[i]).getGolink();
			for (int j = 0; j < clqGOs.size(); j++) {
				if (goitems.contains(strclqGos[j])) {
					count[j]++;
				}
			}
		}

		int max = 0;
		for (int i = 0; i < count.length; i++) {
			// System.out.println(count[i]);
			max = Math.max(count[i], max);
		}

		// System.out.println(max);
		// System.out.println("proteins num:"+pro.length);

		ArrayList<String> mxCommonGoitem = new ArrayList<String>();
		for (int i = 0; i < count.length; i++) {
			if (count[i] == max) {
				mxCommonGoitem.add(strclqGos[i]);
			}
		}
		return mxCommonGoitem;
	}

	private static int DeleteONE(TreeSet<String> clq, ArrayList<Node> proteins,
			HashMap<String, Integer> pro_idx) {

		String[] pros = new String[clq.size()];
		clq.toArray(pros);

		int[] ipros = new int[pros.length];

		for (int i = 0; i < pros.length; i++) {
			ipros[i] = pro_idx.get(pros[i]);
		}

		double orgFitV = CalculateFitness(ipros, proteins);
		// System.out.println("orgFit:"+orgFitV);

		double mxFitV = orgFitV;

		int delIdx = -1;

		for (int i = 0; i < pros.length; i++) {
			double delFitV = CalculateFitness(ipros, proteins, i);
			// System.out.println(delFitV);
			if (delFitV > mxFitV) {
				mxFitV = delFitV;
				delIdx = i;
			}
		}

		// System.out.println("MxFit:"+ mxFitV);
		if (delIdx == -1)
			return -1;
		else
			return ipros[delIdx];
	}

	private static double CalculateFitness(int[] v, ArrayList<Node> proteins) {
		double density = 0;
		double fitV = 0;
		// double[] prefix = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			int[] neighbors2 = proteins.get(v[i]).getNeighbors();
			int degree = neighbors2.length;
			int indegree = intersect(v, neighbors2);
			density += indegree;
			double pi = indegree * 1.0 / degree;
			fitV += pi;
		}
		fitV = fitV * density * 100.0 / v.length / (v.length - 1) / v.length;
		return fitV;
	}

	private static int intersect(int[] a, int[] b) {
		int common = 0;
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < b.length; j++) {
				if (a[i] == b[j])
					common++;
			}
		return common;
	}

	private static double CalculateFitness(int[] v, ArrayList<Node> proteins,
			int n) {
		double density = 0;
		double fitV = 0;
		// double[] prefix = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			if (i == n)
				continue;
			int[] neighbors2 = proteins.get(v[i]).getNeighbors();
			int degree = neighbors2.length;
			int indegree = intersect(v, neighbors2, n);
			density += indegree;
			double pi = indegree * 1.0 / degree;
			fitV += pi;
		}
		fitV = fitV * density * 100.0 / (v.length - 1) / (v.length - 1)
				/ (v.length - 2);
		return fitV;
	}

	private static int intersect(int[] a, int[] b, int n) {
		int common = 0;
		for (int i = 0; i < a.length; i++) {
			if (i == n)
				continue;
			for (int j = 0; j < b.length; j++) {
				if (a[i] == b[j])
					common++;
			}
		}
		return common;
	}

	public static boolean IsBioMeaning(TreeSet<String> clq, double threshold) {
		String[] pro = new String[clq.size()];
		clq.toArray(pro);

		ArrayList<Node> proteins = Network.proteins;
		HashMap<String, Integer> pro_idx = Network.pro_idx;

		HashSet<String> clqGOs = new HashSet<String>();

		int[] ipro = new int[pro.length];
		for (int i = 0; i < pro.length; i++) {
			ipro[i] = pro_idx.get(pro[i]);
			HashSet<String> goitems = proteins.get(ipro[i]).getGolink();
			// System.out.println(pro[i]+":"+goitems.size());
			clqGOs.addAll(goitems);
		}

		int[] count = new int[clqGOs.size()];
		String[] strclqGos = new String[clqGOs.size()];
		clqGOs.toArray(strclqGos);
		// System.out.println("ALL GO ITEMS:"+clqGOs.size());

		for (int i = 0; i < pro.length; i++) {
			HashSet<String> goitems = proteins.get(ipro[i]).getGolink();
			for (int j = 0; j < clqGOs.size(); j++) {
				if (goitems.contains(strclqGos[j])) {
					count[j]++;
				}
			}
		}

		int max = 0;
		for (int i = 0; i < count.length; i++) {
			// System.out.println(count[i]);
			max = Math.max(count[i], max);
		}

		if (max * 1.0 / clq.size() >= threshold)
			return true;
		else
			return false;
	}
	
	private static boolean HasIntersection(TreeSet<Integer> v1, TreeSet<Integer> v2)
	{
		Iterator<Integer> it1 = v1.iterator();
		while(it1.hasNext())
		{
			int i = it1.next();
			Iterator<Integer> it2 = v2.iterator();
			while(it2.hasNext())
			{
				int j= it2.next();
				if(i==j)
					return true;
			}
		}
		return false;
	}
}
