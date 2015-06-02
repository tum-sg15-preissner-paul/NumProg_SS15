import java.net.*;
import java.io.*;
import java.util.*;

public class Crawler {
	private String outFile;
	private String[] startPoints;
	private int maxDepth;
	private String filter;
	private int maxThreads = 10;
	private int threadCount = 0;
	/* Zahl der entdeckten Seiten mit depth<maxDepth */
	private int pageCount = 0;
	/* Zahl der Seiten die bereits gecrawlt wurden */
	private int crawledCount = 0;
	/* enthaelt jede Seite, die gecrawlt oder zum Queue hinzugefuegt wurde */
	private TreeMap<String, Page> pages = new TreeMap<String, Page>();
	private Queue<QueuedPage> q = new LinkedList<QueuedPage>();
	private LinkedList<QueuedPage> currentCrawls = new LinkedList<QueuedPage>();
	private boolean canceled = false, finished = false;

	public Crawler(String outputFile, String[] startURLs, int maxCrawlDepth,
			String URLfilter) {
		outFile = outputFile;
		startPoints = startURLs;
		maxDepth = maxCrawlDepth;
		filter = URLfilter;

		// Queue mit Start-URLs initialisieren
		for (String s : startPoints) {
			s = s.toLowerCase();
			pages.put(s, new Page(s, pageCount));
			pageCount++;
			q.add(new QueuedPage(s, 0));
		}

		startThreads();
	}

	public synchronized void startThreads() {
		QueuedPage qp;

		while (!q.isEmpty() && threadCount < maxThreads && !canceled) {
			qp = q.poll();
			new CrawlerThread(qp, this);
			currentCrawls.add(qp);
			threadCount++;
		}

		// bereite Output vor, falls alle Threads fertig sind und keine weiteren
		// Seite zum Crawlen anstehen
		if ((q.isEmpty() || canceled) && threadCount == 0) {
			LinkMatrix lm = getLinkMatrix();
			try {
				lm.write(outFile);
			} catch (Exception e) {
				System.out.println(outFile
						+ " konnte nicht geschrieben werden.");
			}
			finished = true;
		}
	}

	public void cancel() {
		canceled = true;
	}

	public synchronized void threadFinished(QueuedPage qp,
			LinkedList<String> links) {
		Page currentPage;
		Page temp;
		currentPage = pages.get(qp.url);

		for (String s : links)
			if (s.startsWith(filter)) {
				// System.out.println("[" + qp.depth + "]" + currentPage.url +
				// " -> " + s);
				if (pages.containsKey(s))
					currentPage.links.add(pages.get(s));
				// Es werden nur Seiten beruecksichtigt, die auch gecrawlt
				// werden um PR-Senken zu vermeiden
				else if (qp.depth < maxDepth) {
					temp = new Page(s, pageCount);
					pages.put(s, temp);
					pageCount++;
					currentPage.links.add(temp);
					q.add(new QueuedPage(s, qp.depth + 1));
				}
			}

		currentCrawls.remove(qp);
		threadCount--;
		crawledCount++;
		startThreads();
	}

	public synchronized String getStatus() {
		String s = "";
		if (finished) {
			s += "Crawl abgeschlossen\n";
		} else {
			if (canceled) {
				s += "Es werden keine neuen Seiten gecrawlt. Vorgang wird beendet, "
						+ "wenn alle momentan laufenden Threads fertig sind.\n";
			}
			s += "* Aktive Threads ([Tiefe]URL) *\n";
			for (QueuedPage qp : currentCrawls)
				s += "[" + qp.depth + "]" + qp.url + "\n";
		}
		s += crawledCount + " Seiten gecrawlt";

		return s;
	}

	public boolean hasFinished() {
		return finished;
	}

	private LinkMatrix getLinkMatrix() {
		// erzeuge Matrix
		int L[][] = new int[pageCount][pageCount];
		for (Page j : pages.values()) {
			L[j.nbr][j.nbr] = 1;
			for (Page i : j.links)
				L[i.nbr][j.nbr] = 1;
		}

		// stelle Output zusammen
		LinkMatrix lm = new LinkMatrix();
		lm.L = L;
		lm.urls = new String[pageCount];
		for (Page p : pages.values())
			lm.urls[p.nbr] = p.url;

		return lm;
	}

	private class CrawlerThread implements Runnable {
		private Crawler parent;
		private QueuedPage qp;

		public CrawlerThread(QueuedPage page, Crawler parent) {
			this.parent = parent;
			qp = page;
			(new Thread(this)).start();
		}

		public void run() {
			// Page currentPage;
			LinkedList<String> links;

			// currentPage = pages.get(qp.url);
			try {
				System.out.println(qp.url);

				links = getLinks(readPage(qp.url), qp.url);

				parent.threadFinished(qp, links);
			} catch (Exception e) {
				System.out.println("Exception:" + e.getMessage());
				parent.threadFinished(qp, new LinkedList<String>());
			}
		}
	}

	private static class Page {
		public String url;
		public int nbr;
		public LinkedList<Page> links = new LinkedList<Page>();

		public Page(String u, int n) {
			url = u;
			nbr = n;
		}
	}

	private static class QueuedPage {
		public String url;
		public int depth;

		public QueuedPage(String u, int d) {
			url = u;
			depth = d;
		}
	}

	public static String readPage(String url) throws MalformedURLException,
			IOException {
		int maxSize = 500000;

		URL u = new URL(url);
		URLConnection uc = u.openConnection();
		if (!uc.getContentType().startsWith("text"))
			return "";
		InputStream is = uc.getInputStream();
		BufferedReader br = new BufferedReader(new InputStreamReader(is));

		String s = "";
		String temp = "";
		Boolean firstLine = true;
		while (temp != null && s.length() < maxSize) {
			temp = br.readLine();
			if (temp != null) {
				if (!firstLine)
					s += "\n";
				else
					firstLine = false;
				s += temp;
			}
		}

		return s;
	}

	public static LinkedList<String> getLinks(String data, String url) {
		LinkedList<String> links = new LinkedList<String>();
		URI u;
		String l;
		try {
			u = new URI(url.toLowerCase());
		} catch (Exception e) {
			u = null;
		}

		int pos = 0, posEnd;

		while (true) {
			pos = data.indexOf("<a", pos);
			if (pos == -1)
				break;
			pos = data.indexOf("href", pos + 1);
			if (pos == -1)
				break;
			pos = data.indexOf("\"", pos + 1);
			if (pos == -1)
				break;
			posEnd = data.indexOf("\"", pos + 1);
			if (posEnd == -1)
				break;

			l = data.substring(pos + 1, posEnd).toLowerCase();
			if (u != null)
				l = u.resolve(l).toString();
			links.add(l);
			pos = posEnd + 1;
		}

		return links;
	}
}
