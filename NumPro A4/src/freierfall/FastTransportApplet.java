package freierfall;

import javax.swing.JApplet;

public class FastTransportApplet extends JApplet {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public void init() {
		add(new FastTransportPanel());
		setSize(1150, 600);
		setVisible(true);
	}
}
