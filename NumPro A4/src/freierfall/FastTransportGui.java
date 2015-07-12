package freierfall;

import javax.swing.JFrame;

public class FastTransportGui extends JFrame {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public FastTransportGui() {
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		add(new FastTransportPanel());
		pack();
	}

	public static void main(String[] args) {
		FastTransportGui frame = new FastTransportGui();
		frame.setVisible(true);
	}
}
