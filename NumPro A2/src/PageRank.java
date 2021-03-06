import java.util.Arrays;
import java.util.Comparator;

public class PageRank {

	/**
	 * Diese Methode erstellt die Matrix A~ fuer das PageRank-Verfahren
	 * PARAMETER: 
	 * L: die Linkmatrix (s. Aufgabenblatt) 
	 * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
	 *      zufaellig irgendeine Seite zu besuchen
	 */
	public static double[][] buildProbabilityMatrix(int[][] L, double rho) {
		int n = L.length;
		double[][] A = new double[n][n];
		double[] linkCounts = new double[n];
		
		//Lij = 1 if j links to i
		//for each node (j), calculate how many outgoing links
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				 if(L[i][j] == 1) linkCounts[j]++;
			}
		}
		
		//-> Aix = 1/(sum of all Ljx), if j links to i (!!!)
		//~aij = (1-p) * aij + p/n ---- p = rho
		//Aij = ((1-p) * (1/sum)) + rho/n;
		double nD = (double) n;
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				//aij is only 1/sum if j links to i, so if Lij == 1
				if(L[i][j] == 1) {
					A[i][j] = ((1 - rho) * (1.0/linkCounts[j]));
				}
				//we still add p/n, since that isnt affected by aij
				A[i][j] += rho/nD;
			}
		}
		
		return A;
	}

	/**
	 * Diese Methode berechnet die PageRanks der einzelnen Seiten,
	 * also das Gleichgewicht der Aufenthaltswahrscheinlichkeiten.
	 * (Entspricht dem p-Strich aus der Angabe)
	 * Die Ausgabe muss dazu noch normiert sein.
	 * PARAMETER:
	 * L: die Linkmatrix (s. Aufgabenblatt) 
	 * rho: Wahrscheinlichkeit, zufaellig irgendeine Seite zu besuchen
	 * ,anstatt einem Link zu folgen.
	 *      
	 */
	public static double[] rank(int[][] L, double rho) {
		int n = L.length;
		double[][] probMatrix = buildProbabilityMatrix(L, rho);
		double[] p = new double[n];
		
		//since we want the solution of (~A - I)p = 0, deduct I from ~A
		for(int i = 0; i < n; i++) {
			probMatrix[i][i] -= 1.0;
		}
		
		p = Gauss.solveSing(probMatrix);
		
		//normalize p
		double lamba = 0;
		for(int i = 0; i < n; i++) {
			lamba += p[i];
		}
		for(int i = 0; i < n; i++) {
			p[i] = p[i] / lamba;
		}
		
		return p;
	}

	/**
	 * Diese Methode erstellt eine Rangliste der uebergebenen URLs nach
	 * absteigendem PageRank. 
 	 * PARAMETER:
 	 * urls: Die URLs der betrachteten Seiten
 	 * L: die Linkmatrix (s. Aufgabenblatt) 
 	 * rho: Wahrscheinlichkeit, anstatt einem Link zu folgen,
 	 *      zufaellig irgendeine Seite zu besuchen
	 */ 
	public static String[] getSortedURLs(String[] urls, int[][] L, double rho) {
		int n = L.length;

		double[] p = rank(L, rho);

		RankPair[] sortedPairs = new RankPair[n];
		for (int i = 0; i < n; i++) {
			sortedPairs[i] = new RankPair(urls[i], p[i]);
		}

		Arrays.sort(sortedPairs, new Comparator<RankPair>() {

			@Override
			public int compare(RankPair o1, RankPair o2) {
				return -Double.compare(o1.pr, o2.pr);
			}
		});

		String[] sortedUrls = new String[n];
		for (int i = 0; i < n; i++) {
			sortedUrls[i] = sortedPairs[i].url;
		}

		return sortedUrls;
	}

	/**
	 * Ein RankPair besteht aus einer URL und dem zugehoerigen Rang, und dient
	 * als Hilfsklasse zum Sortieren der Urls
	 */
	private static class RankPair {
		public String url;
		public double pr;

		public RankPair(String u, double p) {
			url = u;
			pr = p;
		}
	}
}
