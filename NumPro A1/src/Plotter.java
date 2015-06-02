import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import javax.swing.JPanel;

/**
 * @author Christoph Riesinger (riesinge@in.tum.de)
 * @since November 06, 2011
 * @version 1.0
 * 
 *          This class is a simple implementation of a plotter. It accepts two
 *          arrays of float values which represent the x- and y-coordinates of
 *          points which should be plotted in a carthesian coordinate system. x-
 *          and y-axis are scaled logarithmically. No interpolation between the
 *          dots is done!
 */
public class Plotter extends JPanel {
	/* Just to avoid compiler warning. */
	private static final long serialVersionUID = -8725968470735352529L;

	/* Top, bottom, left and right margin of the coordinate system in the plot. */
	private final int PADDING = 20;

	/* Class members. Are set in class constructor. */
	private float[] xData = null;
	private float[] yData = null;
	private float minX = Float.MAX_VALUE;
	private float minY = Float.MAX_VALUE;
	private float maxX = -Float.MAX_VALUE;
	private float maxY = -Float.MAX_VALUE;

	/**
	 * Constructor of this class. Assigns the passed x- and y-values of the
	 * points to plot to the internal private member variables.
	 * 
	 * @param xData
	 *            x-values of the points which should be plotted by this class.
	 * @param yData
	 *            y-values of the points which should be plotted by this class.
	 * @throws InstantiationException
	 *             The lengths of the x- and y-coorinates arrays have to be
	 *             equal. Elsewise an exception is thrown.
	 */
	public Plotter(float[] xData, float[] yData) throws InstantiationException {
		/*
		 * Make sure the arrays which contain the x- and y-ccordinates which
		 * should be plotted by this class have the same length.
		 */
		if (xData.length != yData.length) {
			throw (new InstantiationException(
					"The arrays for the x- and y-components of the "
							+ "coordinates have to be of the same length."));
		}

		this.xData = xData;
		this.yData = yData;

		/*
		 * Determine the smallest and largest value which should be plotted by
		 * this class. These values are the boundaries of the axes of the
		 * coordinate system which will be plotted.
		 */
		for (int i = 0; i < xData.length; i++) {
			if (xData[i] < minX) {
				minX = xData[i];
			}
			if (xData[i] > maxX) {
				maxX = xData[i];
			}
			if (yData[i] < minY) {
				minY = yData[i];
			}
			if (yData[i] > maxY) {
				maxY = yData[i];
			}
		}
		if (1.0d / Math.sqrt(maxX) < minY) {
			minY = (float) (1.0d / Math.sqrt(maxX));
		}
		if (1.0d / Math.sqrt(minX) < minY) {
			minY = (float) (1.0d / Math.sqrt(minX));
		}
		if (1.0d / Math.sqrt(maxX) > maxY) {
			maxY = (float) (1.0d / Math.sqrt(maxX));
		}
		if (1.0d / Math.sqrt(minX) > maxY) {
			maxY = (float) (1.0d / Math.sqrt(minX));
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see javax.swing.JComponent#paintComponent(java.awt.Graphics)
	 */
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);

		Graphics2D graphics = (Graphics2D) g;
		graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		int height = getHeight();
		int width = getWidth();
		float x, y;

		/* draw x-axis */
		graphics.draw(new Line2D.Double(PADDING, PADDING, PADDING, height
				- PADDING));
		graphics.draw(new Line2D.Double(
				PADDING + 0.25f * (width - 2 * PADDING), height - PADDING - 5,
				PADDING + 0.25f * (width - 2 * PADDING), height - PADDING + 5));
		graphics.draw(new Line2D.Double(0.5f * width, height - PADDING - 5,
				0.5f * width, height - PADDING + 5));
		graphics.draw(new Line2D.Double(
				PADDING + 0.75f * (width - 2 * PADDING), height - PADDING - 5,
				PADDING + 0.75f * (width - 2 * PADDING), height - PADDING + 5));
		graphics.draw(new Line2D.Double(width - PADDING, height - PADDING - 5,
				width - PADDING, height - PADDING + 5));

		/* draw x-axis caption */
		graphics.drawString((new Float(minX)).toString(), PADDING + 2, height
				- PADDING + 12);
		graphics.drawString(
				(new Float(Math.pow(10.0d, 0.25f * Math.log10(maxX) + 0.75f
						* Math.log10(minX)))).toString(), PADDING + 0.25f
						* (width - 2 * PADDING) + 2, height - PADDING + 12);
		graphics.drawString(
				(new Float(Math.pow(10.0d,
						0.5f * (Math.log10(maxX) + Math.log10(minX)))))
						.toString(), 0.5f * width + 2, height - PADDING + 12);
		graphics.drawString(
				(new Float(Math.pow(10.0d, 0.75f * Math.log10(maxX) + 0.25f
						* Math.log10(minX)))).toString(), PADDING + 0.75f
						* (width - 2 * PADDING) + 2, height - PADDING + 12);
		graphics.drawString((new Float(maxX)).toString(), width - PADDING + 2,
				height - PADDING + 12);

		/* draw y-axis */
		graphics.draw(new Line2D.Double(PADDING, height - PADDING, width
				- PADDING, height - PADDING));
		graphics.draw(new Line2D.Double(PADDING - 5, height
				- (PADDING + 0.25f * (height - 2 * PADDING)), PADDING + 5,
				height - (PADDING + 0.25f * (height - 2 * PADDING))));
		graphics.draw(new Line2D.Double(PADDING - 5, 0.5f * height,
				PADDING + 5, 0.5f * height));
		graphics.draw(new Line2D.Double(PADDING - 5, height
				- (PADDING + 0.75f * (height - 2 * PADDING)), PADDING + 5,
				height - (PADDING + 0.75f * (height - 2 * PADDING))));
		graphics.draw(new Line2D.Double(PADDING - 5, height
				- (PADDING + (height - 2 * PADDING)), PADDING + 5, height
				- (PADDING + (height - 2 * PADDING))));

		/* draw y-axis caption */
		graphics.drawString((new Float(minY)).toString(), PADDING + 2,
				(height - 2 * PADDING) + PADDING - 2);
		graphics.drawString(
				(new Float(Math.pow(10.0d, 0.25f * Math.log10(maxY) + 0.75f
						* Math.log10(minY)))).toString(), PADDING + 2, height
						- (PADDING + 0.25f * (height - 2 * PADDING)) - 2);
		graphics.drawString(
				(new Float(Math.pow(10.0d,
						0.5f * (Math.log10(maxY) + Math.log10(minY)))))
						.toString(), PADDING + 2, 0.5f * height - 2);
		graphics.drawString(
				(new Float(Math.pow(10.0d, 0.75f * Math.log10(maxY) + 0.25f
						* Math.log10(minY)))).toString(), PADDING + 2, height
						- (PADDING + 0.75f * (height - 2 * PADDING)) - 2);
		graphics.drawString((new Float(maxY)).toString(), PADDING + 2,
				PADDING - 2);

		/* draw "exact" solution */
		graphics.setPaint(Color.GREEN);

		for (int i = 0; i < xData.length; i++) {
			x = scaleX(xData[i], width);
			y = scaleY((float) (1.0d / Math.sqrt(xData[i])), height);
			graphics.fill(new Ellipse2D.Float(x - 1.0f, y - 1.0f, 2.0f, 2.0f));
		}

		/* draw assigned values */
		graphics.setPaint(Color.red);

		for (int i = 0; i < xData.length; i++) {
			x = scaleX(xData[i], width);
			y = scaleY(yData[i], height);
			graphics.fill(new Ellipse2D.Float(x - 1.0f, y - 1.0f, 2.0f, 2.0f));
		}
	}

	/**
	 * The x-values are logarithmically scaled before they are drawn. This is
	 * done by this method. The padding of the coordinate system is respected.
	 * 
	 * @param x
	 *            Value which should be scaled logarithmically.
	 * @param width
	 *            Width of the plottable area.
	 * @return Logarithmically scaled value of x.
	 */
	private float scaleX(float x, float width) {
		float xScale = (float) (width - 2 * PADDING)
				/ (float) (Math.log10(maxX) - Math.log10(minX));
		float xOffset = (float) -Math.log10(minX);
		float result = PADDING + xScale * (xOffset + (float) Math.log10(x));

		return result;
	}

	/**
	 * The y-values are logarithmically scaled before they are drawn. This is
	 * done by this method. The padding of the coordinate system is respected.
	 * 
	 * @param y
	 *            Value which should be scaled logarithmically.
	 * @param height
	 *            Height of the plottable area.
	 * @return Logarithmically scaled value of y.
	 */
	private float scaleY(float y, float height) {
		float yScale = (float) (height - 2 * PADDING)
				/ (float) (Math.log10(maxY) - Math.log10(minY));
		float yOffset = (float) -Math.log10(minY);
		float result = height
				- (PADDING + yScale * (yOffset + (float) Math.log10(y)));

		return result;
	}
}
