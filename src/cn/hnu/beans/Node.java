package cn.hnu.beans;

import java.util.HashSet;

public class Node {
	
	private int idx;
	private String name;
	private HashSet<String> golink;
	private int[] neighbors;
	private int degree;
	
	public Node(String name, int idx){
		this.setName(name);
		this.setIdx(idx);
		golink = new HashSet<String>();
		degree = 0;
	}
	
	public void setGolink(HashSet<String> golink) {
		this.golink = golink;
	}
	public HashSet<String> getGolink() {
		return golink;
	}

	public void setIdx(int idx) {
		this.idx = idx;
	}

	public int getIdx() {
		return idx;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}

	public void setDegree(int degree) {
		this.degree = degree;
	}

	public int getDegree() {
		return degree;
	}

	public void setNeighbors(int[] neighbors) {
		this.neighbors = neighbors;
	}

	public int[] getNeighbors() {
		return neighbors;
	}
	
}
