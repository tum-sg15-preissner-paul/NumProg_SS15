import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JMenuBar;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JTextArea;

import java.awt.GridBagLayout;
import javax.swing.JTextPane;
import java.awt.GridBagConstraints;
import java.util.*;
import javax.swing.JButton;
import javax.swing.JLabel;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import javax.swing.SwingConstants;
import javax.swing.JTextField;
import java.awt.Font;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.JScrollPane;
import javax.swing.filechooser.FileFilter;

public class GUI extends JFrame {
	private Crawler crawler;
	private Timer timer;

	private static final long serialVersionUID = 1L;

	private JPanel jContentPane = null;

	private JTextPane jtpCrawls = null;

	private JButton jbGo = null;
	private JButton jbCancel = null;
	private JLabel jlStartURLs = null;
	private JTextPane jtpStartURLs = null;
	private JLabel jlFilter = null;
	private JTextField jtfFilter = null;
	private JLabel jlFilterExp = null;
	private JLabel jlDepth = null;
	private JTextField jtfMaxDepth = null;
	private JLabel jlOutput = null;
	private JTextField jtfOutput = null;
	private JScrollPane jspStartURLs = null;
	private JScrollPane jspCrawls = null;

	/**
	 * This is the default constructor
	 */
	public GUI() {
		super();
		initialize();
	}

	/**
	 * This method initializes this
	 * 
	 * @return void
	 */
	private void initialize() {
		this.setSize(414, 431);
		this.setLocation(86, 200);
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setContentPane(getJContentPane());
		this.setTitle("Crawler");
	}

	/**
	 * This method initializes jContentPane
	 * 
	 * @return javax.swing.JPanel
	 */
	private JPanel getJContentPane() {
		if (jContentPane == null) {
			GridBagConstraints gridBagConstraints13 = new GridBagConstraints();
			gridBagConstraints13.fill = GridBagConstraints.BOTH;
			gridBagConstraints13.gridy = 8;
			gridBagConstraints13.weightx = 1.0;
			gridBagConstraints13.weighty = 4.0D;
			gridBagConstraints13.gridwidth = 2;
			gridBagConstraints13.gridx = 0;
			GridBagConstraints gridBagConstraints12 = new GridBagConstraints();
			gridBagConstraints12.fill = GridBagConstraints.BOTH;
			gridBagConstraints12.gridy = 1;
			gridBagConstraints12.weightx = 1.0;
			gridBagConstraints12.weighty = 1.0;
			gridBagConstraints12.gridwidth = 3;
			gridBagConstraints12.gridx = 0;
			GridBagConstraints gridBagConstraints111 = new GridBagConstraints();
			gridBagConstraints111.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints111.gridy = 6;
			gridBagConstraints111.weightx = 1.0;
			gridBagConstraints111.gridx = 1;
			GridBagConstraints gridBagConstraints10 = new GridBagConstraints();
			gridBagConstraints10.gridx = 0;
			gridBagConstraints10.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints10.gridy = 6;
			jlOutput = new JLabel();
			jlOutput.setText("Ausgabedatei: ");
			jlOutput.setHorizontalAlignment(SwingConstants.RIGHT);
			GridBagConstraints gridBagConstraints9 = new GridBagConstraints();
			gridBagConstraints9.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints9.gridy = 5;
			gridBagConstraints9.weightx = 1.0;
			gridBagConstraints9.gridx = 1;
			GridBagConstraints gridBagConstraints7 = new GridBagConstraints();
			gridBagConstraints7.gridx = 0;
			gridBagConstraints7.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints7.gridy = 5;
			jlDepth = new JLabel();
			jlDepth.setText("max. Suchtiefe: ");
			jlDepth.setHorizontalTextPosition(SwingConstants.RIGHT);
			jlDepth.setHorizontalAlignment(SwingConstants.RIGHT);
			GridBagConstraints gridBagConstraints6 = new GridBagConstraints();
			gridBagConstraints6.gridx = 0;
			gridBagConstraints6.gridwidth = 2;
			gridBagConstraints6.gridy = 4;
			jlFilterExp = new JLabel();
			jlFilterExp
					.setText("(Es werden nur URLs gecrawlt, die mit FILTER beginnen)");
			jlFilterExp.setFont(new Font("Dialog", Font.PLAIN, 12));
			GridBagConstraints gridBagConstraints5 = new GridBagConstraints();
			gridBagConstraints5.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints5.gridy = 3;
			gridBagConstraints5.weightx = 0.0D;
			gridBagConstraints5.gridwidth = 1;
			gridBagConstraints5.gridx = 1;
			GridBagConstraints gridBagConstraints4 = new GridBagConstraints();
			gridBagConstraints4.gridx = 0;
			gridBagConstraints4.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints4.gridwidth = 1;
			gridBagConstraints4.weightx = 0.0D;
			gridBagConstraints4.gridy = 3;
			jlFilter = new JLabel();
			jlFilter.setText(" Filter: ");
			jlFilter.setHorizontalAlignment(SwingConstants.RIGHT);
			jlFilter.setHorizontalTextPosition(SwingConstants.RIGHT);
			GridBagConstraints gridBagConstraints2 = new GridBagConstraints();
			gridBagConstraints2.gridx = 0;
			gridBagConstraints2.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints2.gridwidth = 2;
			gridBagConstraints2.gridy = 0;
			jlStartURLs = new JLabel();
			jlStartURLs.setText(" StartURLs (eine URL pro Zeile)");
			jlStartURLs.setHorizontalAlignment(SwingConstants.LEFT);
			jlStartURLs.setDisplayedMnemonic(KeyEvent.VK_UNDEFINED);
			GridBagConstraints gridBagConstraints11 = new GridBagConstraints();
			gridBagConstraints11.gridx = 1;
			gridBagConstraints11.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints11.weightx = 1.0D;
			gridBagConstraints11.gridy = 7;
			GridBagConstraints gridBagConstraints1 = new GridBagConstraints();
			gridBagConstraints1.gridx = 0;
			gridBagConstraints1.fill = GridBagConstraints.HORIZONTAL;
			gridBagConstraints1.weightx = 1.0D;
			gridBagConstraints1.gridwidth = 1;
			gridBagConstraints1.gridy = 7;
			jContentPane = new JPanel();
			jContentPane.setLayout(new GridBagLayout());
			jContentPane.add(getJbGo(), gridBagConstraints1);
			jContentPane.add(getJbCancel(), gridBagConstraints11);
			jContentPane.add(jlStartURLs, gridBagConstraints2);
			jContentPane.add(jlFilter, gridBagConstraints4);
			jContentPane.add(getJtfFilter(), gridBagConstraints5);
			jContentPane.add(jlFilterExp, gridBagConstraints6);
			jContentPane.add(jlDepth, gridBagConstraints7);
			jContentPane.add(getJtfMaxDepth(), gridBagConstraints9);
			jContentPane.add(jlOutput, gridBagConstraints10);
			jContentPane.add(getJtfOutput(), gridBagConstraints111);
			jContentPane.add(getJspStartURLs(), gridBagConstraints12);
			jContentPane.add(getJspCrawls(), gridBagConstraints13);
		}
		return jContentPane;
	}

	/**
	 * This method initializes jtpCrawls
	 * 
	 * @return javax.swing.JTextPane
	 */
	private JTextPane getJtpCrawls() {
		if (jtpCrawls == null) {
			jtpCrawls = new JTextPane();
			jtpCrawls.setEditable(false);
		}
		return jtpCrawls;
	}

	class updateStatus extends TimerTask {
		// private int i = 0;

		public void run() {
			if (crawler != null) {
				jtpCrawls.setText(crawler.getStatus());
				if (crawler.hasFinished()) {
					timer.cancel();
					jbGo.setEnabled(true);
				}
				// System.out.println(i); i++;
			}
		}

	}

	/**
	 * This method initializes jbGo
	 * 
	 * @return javax.swing.JButton
	 */
	private JButton getJbGo() {
		if (jbGo == null) {
			jbGo = new JButton();
			jbGo.setText("Los!");
			jbGo.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					crawler = new Crawler(jtfOutput.getText(), jtpStartURLs
							.getText().split("\n"), Integer.valueOf(jtfMaxDepth
							.getText()), jtfFilter.getText());
					timer = new Timer();
					timer.schedule(new updateStatus(), 1000, 1000);
					jbGo.setEnabled(false);
					jbCancel.setEnabled(true);
				}
			});
		}
		return jbGo;
	}

	/**
	 * This method initializes jbCancel
	 * 
	 * @return javax.swing.JButton
	 */
	private JButton getJbCancel() {
		if (jbCancel == null) {
			jbCancel = new JButton();
			jbCancel.setText("Abbrechen");
			jbCancel.setEnabled(false);
			jbCancel.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					crawler.cancel();
					jbCancel.setEnabled(false);
				}
			});
		}
		return jbCancel;
	}

	/**
	 * This method initializes jtpStartURLs
	 * 
	 * @return javax.swing.JTextPane
	 */
	private JTextPane getJtpStartURLs() {
		if (jtpStartURLs == null) {
			jtpStartURLs = new JTextPane();
			jtpStartURLs
					.setText("http://www-e.uni-magdeburg.de/harbich/webcrawling.php");
		}
		return jtpStartURLs;
	}

	/**
	 * This method initializes jtfFilter
	 * 
	 * @return javax.swing.JTextField
	 */
	private JTextField getJtfFilter() {
		if (jtfFilter == null) {
			jtfFilter = new JTextField();
			jtfFilter.setText("http://");
		}
		return jtfFilter;
	}

	/**
	 * This method initializes jtfMaxDepth
	 * 
	 * @return javax.swing.JTextField
	 */
	private JTextField getJtfMaxDepth() {
		if (jtfMaxDepth == null) {
			jtfMaxDepth = new JTextField();
			jtfMaxDepth.setHorizontalAlignment(JTextField.RIGHT);
			jtfMaxDepth.setText("3");
		}
		return jtfMaxDepth;
	}

	/**
	 * This method initializes jtfOutput
	 * 
	 * @return javax.swing.JTextField
	 */
	private JTextField getJtfOutput() {
		if (jtfOutput == null) {
			jtfOutput = new JTextField();
			jtfOutput.setText("output.txt");
		}
		return jtfOutput;
	}

	/**
	 * This method initializes jspStartURLs
	 * 
	 * @return javax.swing.JScrollPane
	 */
	private JScrollPane getJspStartURLs() {
		if (jspStartURLs == null) {
			jspStartURLs = new JScrollPane();
			jspStartURLs.setViewportView(getJtpStartURLs());
		}
		return jspStartURLs;
	}

	/**
	 * This method initializes jspCrawls
	 * 
	 * @return javax.swing.JScrollPane
	 */
	private JScrollPane getJspCrawls() {
		if (jspCrawls == null) {
			jspCrawls = new JScrollPane();
			jspCrawls.setViewportView(getJtpCrawls());
		}
		return jspCrawls;
	}

	public static void main(String[] args) {
		(new GUI()).setVisible(true);
		(new Search()).setVisible(true);
	}

} // @jve:decl-index=0:visual-constraint="10,10"

class Search extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;

	JMenuBar mbar = new JMenuBar();

	JButton openbutton = new JButton("Linkmatrix laden");
	JLabel rholabel = new JLabel("Rho:");
	JTextField rhofield = new JTextField("0.15");

	JTextArea textarea = new JTextArea(
			"Hier steht die Rank-Liste der Webseiten", 20, 50);

	public Search() {
		super("Search");

		mbar.add(openbutton);
		mbar.add(rholabel);
		mbar.add(rhofield);

		openbutton.addActionListener(this);

		setJMenuBar(mbar);
		add(textarea);

		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(500, 500);
		setLocation(500, 200);
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		JComponent source = (JComponent) (e.getSource());

		if (source == (JComponent) openbutton) {
			JFileChooser d = new JFileChooser();
			d.setCurrentDirectory(new File("./"));

			d.setFileFilter(new TxtFilter());
			int returnVal = d.showOpenDialog(null);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				try {

					File file = d.getSelectedFile();
					LinkMatrix lm = new LinkMatrix();
					lm.read(file.getAbsolutePath());
					double rho = Double.parseDouble(rhofield.getText());
					double[] rank = PageRank.rank(lm.L, rho);
					Arrays.sort(rank);
					String r[] = PageRank.getSortedURLs(lm.urls, lm.L, rho);
					textarea.setText("");
					for (int i = 0; i < r.length && i < 20; i++) {
						textarea.append((int) (10000 * rank[rank.length - i - 1])
								/ 100.0 + "%\t");
						textarea.append(r[i] + "\n");
					}

					System.out.println("Opening: " + file.getName() + ".\n");
				} catch (FileNotFoundException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}

			} else {
				System.out.println("Open command cancelled by user.\n");
			}
		} else if (source == (JComponent) rhofield) {

		}

	}

	class TxtFilter extends FileFilter {

		@Override
		public boolean accept(File f) {
			return f.isDirectory()
					|| f.getName().toLowerCase().endsWith(".txt");
		}

		@Override
		public String getDescription() {
			return "*.txt";
		}

	}

}