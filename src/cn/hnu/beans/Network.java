package cn.hnu.beans;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

import cn.hnu.actions.FileOperation;
import cn.hnu.actions.GoAnnotation;
import cn.hnu.actions.PPNHandler;



public class Network {
	
	public static int ProteinNUM ; //蛋白质网络中的蛋白质数目 4930
	
	public static int PPIsNUM;
	
	public static ArrayList<Node> proteins;  //网络中的蛋白质
	
//	public static HashMap<Integer,Node> idx_node;  ////
	public static HashMap<String,Integer> pro_idx; //蛋白质到编号id
	
//	public static HashMap<String,Node> pro_node;
	
	public static int graph[][]; //图邻接矩阵
	
	public static HashMap<String, HashSet<Integer>> goitem_pros;//go items对应的蛋白质集合
	
	
	public static int GoitemsNUM ; // GO items条数
	
	public static double[][] goitemssim; // GO items之间的功能相似性
	
	public static HashMap<String, Integer> goitems_idx; //GO items映射对应的id
	
	
	public static double pro_gosim[][]; //蛋白质之间的GO相似性
	
}
