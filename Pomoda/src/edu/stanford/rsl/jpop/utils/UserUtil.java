package edu.stanford.rsl.jpop.utils;

import javax.swing.JOptionPane;


/**
 * Class for obtaining information from the user easily.
 * 
 * @author akmaier
 *
 */
public abstract class UserUtil {
	
	/**
	 * Queries the User for an Integer value using Swing.
	 * @param message
	 * @param initialValue
	 * @return
	 * @throws Exception
	 */
	public static int queryInt(String message, int initialValue) throws Exception{
		String input = JOptionPane.showInputDialog(message, "" + initialValue);
		if (input == null) throw new Exception("Selection aborted");
		return Integer.parseInt(input);
	}
	
	/**
	 * Queries the User for a Double values using Swing.
	 * @param message
	 * @param initialValue
	 * @return
	 * @throws Exception
	 */
	public static double queryDouble(String message, double initialValue) throws Exception{
		String input = JOptionPane.showInputDialog(message, "" + initialValue);
		if (input == null) throw new Exception("Selection aborted");
		return Double.parseDouble(input);
	}
	
	/**
	 * Queries the User for a String value.
	 * @param message
	 * @param initialValue
	 * @return
	 * @throws Exception
	 */
	public static String queryString(String message, String initialValue) throws Exception{
		String input = JOptionPane.showInputDialog(message, "" + initialValue);
		if (input == null) throw new Exception("Selection aborted");
		return input;
	}
	

	/**
	 * Queries the User for a Boolean value.
	 * @param message
	 * @return
	 * @throws Exception
	 */
	public static boolean queryBoolean(String message) throws Exception{
		int revan = JOptionPane.showConfirmDialog(null, message);
		return (revan == JOptionPane.YES_OPTION);
	}
	
	/**
	 * Asks the User to select an Object from a given array of Objects.
	 * @param message
	 * @param messageTitle
	 * @param objects
	 * @param initialObject
	 * @return
	 * @throws Exception
	 */
	public static Object chooseObject(String message, String messageTitle, Object [] objects, Object initialObject) throws Exception{
		Object input = JOptionPane.showInputDialog(null, message, messageTitle, JOptionPane.INFORMATION_MESSAGE, null, objects, initialObject);
		if (input == null) throw new Exception("Selection aborted");
		return input;
	}
	
}
