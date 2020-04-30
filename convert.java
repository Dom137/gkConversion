package com.OpenAdvice.MTA;

import java.util.ArrayList;
import java.util.Properties;
import org.apache.log4j.Logger;

/**
* Klasse zur Durchfuehrung von Koordinatentransformationen
* @seeAlso Map projections: a reference manual Von Lev M. Bugayevskiy,John Parr Snyder
*     		  http://books.google.de/books?id=vTLAqGTAc8cC&lpg=PP1&ots=tHB3YyB59B&dq=map%20projections%20reference%20manual&pg=PP1#v=onepage&q=&f=false
* @ImplementierteTransformationen
* Umrechung von Gauss-Krüger nach WGS84 nach Helmert (7-Parameter-Transformation)
*    @seeAlso http://de.wikipedia.org/wiki/Helmert-Transformation
*/
public class ConvertCoords {
	private Properties m_properties;
	private static Logger log = null;
	ArrayList<Double> returnArray = new ArrayList<Double>();
	
	/***************************************************************************************************************/
	/**
	 * @param props: Properties Map, in der sich mindestens angaben über
	 *  - den Namen eines zu verwendenten Loggers befinden
	 *  - alle nötigen Parameter für eine Transformation zu finden sind
	 */
	public ConvertCoords(Properties props){		
		m_properties = props;
		log = Logger.getLogger(m_properties.getProperty("LoggerName"));	
	}
	
	/***************************************************************************************************************/
	/**
	 * Umrechung von Gauss-Krüger nach WGS84 nach Helmert (7-Parameter-Transformation)
	 * @param Rechtswert
	 * @param Hochwert
	 * @return Arraylist: [0] = longitude
	 * 					  [1] = latitude
	 */	
	public ArrayList transformGK2WGS84(double Rechtswert, double Hochwert) {
		log.debug("transformGK2WGS84():Funktionsaufruf");
		log.debug("transformGK2WGS84():Rechtswert: "+Rechtswert+" Hochwert: "+Hochwert);
		returnArray = new ArrayList<Double>();
		//+++++++++++++++++++++++++++++++++++++++++Init+++++++++++++++++++++++++++++++++++++++++++++++
		double Latitude;
		double Longitude;
		double Altitude;

		double hochwert = Hochwert;
		double rechtswert = Rechtswert;
		double h = 0;
		//Zone finden: Ist immer die erste Stelle des Rechtswertes
		String a = String.valueOf(rechtswert).substring(0, 1);
        int zone = Integer.valueOf(a);
        int bMeridian = zone * 3;  		
        
  
		String helmertGebiet = "Deutschland Potsdamm 2001";
		// Ellepsoid-Koordinaten auf dem Bessel Ellipsoid
		double B; //elliptische Breite
		double L; //elliptische laenge
		// Vektoren in DHDN/Bessel = Deutsche Hauptdreiecksnetz
		// zur Berechnung des Referenzellipsoides
		// http://de.wikipedia.org/wiki/Referenzellipsoid
		double xB;
		double yB;
		double zB;
		// Vektoren in WGS84
		double xW;
		double yW;
		double zW;
		//Konstanten
		//WGS84 Ellipsoid
        double aW = Double.valueOf(m_properties.getProperty("WGS84_gr_Halbachse")); //6378137=große Halbachse
        double bW = Double.valueOf(m_properties.getProperty("WGS84_kl_Halbachse")); //6356752.3141 = kleine Halbachse
        double e2W = (Math.pow(aW, 2) - Math.pow(bW, 2)) / Math.pow(aW, 2); //1.Numerische Exzentrität
        //Bessel Ellipsoid
        //http://de.wikipedia.org/wiki/Bessel-Ellipsoid
        double aB = Double.valueOf(m_properties.getProperty("BesselAequatorAchse")); // 6377397.155 = aequatorial Achse
        double bB = Double.valueOf(m_properties.getProperty("BesselPolarAchse"));// 6356078.962 = polar Achse
        double e2B = (aB * aB - bB * bB) / (aB * aB);		
		// Helmert-Parameter
        // Deutschland von WGS84 nach DNDH/Potsdamm2001    
        // Parameter siehe Wikipedia: http://de.wikipedia.org/wiki/Helmert-Transformation
        // Achtung: inverse Transformation => bei allen Parametern Vorzeichen umdrehen
		double dx = Double.valueOf(m_properties.getProperty("HelmertTransX")); //598.1 = Translation in X
		double dy = Double.valueOf(m_properties.getProperty("HelmertTransY")); // 73.7 = Translation in Y 
		double dz = Double.valueOf(m_properties.getProperty("HelmertTransZ")); //418.2 = Translation in Z
		double ex = Double.valueOf( m_properties.getProperty("HelmertDrehwinkelX")); //-0.202 = Drehwinkel in Bogensekunden un die x-Achse
		double ey = Double.valueOf( m_properties.getProperty("HelmertDrehwinkelY")); //-0.045 = Drehwinkel in Bogensekunden un die y-Achse
		double ez = Double.valueOf( m_properties.getProperty("HelmertDrehwinkelZ")); //2.455 = Drehwinkel in Bogensekunden un die z-Achse
		double m = Double.valueOf( m_properties.getProperty("HelmertMassstabPPM")); //6.7 = Maßstabsfaktor in ppm ;
		//+++++++++++++++++++++++++++++++++++++++++Init Ende+++++++++++++++++++++++++++++++++++++++++++
		
		//Bessel Ellipsoid
        double n = (aB - bB) / (aB + bB);
        double alpha = (aB + bB) / 2.0 * (1.0 + 1.0 / 4.0 * n * n + 1.0 / 64.0 * Math.pow(n, 4));
        double beta = 3.0 / 2.0 * n - 27.0 / 32.0 * Math.pow(n, 3) + 269.0 / 512.0 * Math.pow(n, 5);
        double gamma = 21.0 / 16.0 * n * n - 55.0 / 32.0 * Math.pow(n, 4);
        double delta = 151.0 / 96.0 * Math.pow(n, 3) - 417.0 / 128.0 * Math.pow(n, 5);
        double epsilon = 1097.0 / 512.0 * Math.pow(n, 4);
        double y0 = bMeridian / 3.0;
        double y = rechtswert - y0 * 1000000 - 500000;
        double B0 = hochwert / alpha;
        double Bf = B0 + beta * Math.sin(2 * B0) + gamma * Math.sin(4 * B0) + delta * Math.sin(6 * B0) + epsilon * Math.sin(8 * B0);
        double Nf = aB / Math.sqrt(1.0 - e2B * Math.pow(Math.sin(Bf), 2));
        double ETAf = Math.sqrt((aB * aB) / (bB * bB) * e2B * Math.pow(Math.cos(Bf), 2));
        double tf = Math.tan(Bf);
        double b1 = tf / 2.0 / (Nf * Nf) * (-1.0 - (ETAf * ETAf)) * (y * y);
        double b2 = tf / 24.0 / Math.pow(Nf, 4) * (5.0 + 3.0 * (tf * tf) + 6.0 * (ETAf * ETAf) - 6.0 * (tf * tf) * (ETAf * ETAf) - 4.0 * Math.pow(ETAf, 4) - 9.0 * (tf * tf) * Math.pow(ETAf, 4)) * Math.pow(y, 4);
        B = (Bf + b1 + b2) * 180 / Math.PI;
        double l1 = 1.0 / Nf / Math.cos(Bf) * y;
        double l2 = 1.0 / 6.0 / Math.pow(Nf, 3) / Math.cos(Bf) * (-1.0 - 2.0 * (tf * tf) - (ETAf * ETAf)) * Math.pow(y, 3);
        L = bMeridian + (l1 + l2) * 180 / Math.PI;
        
        //Ellipsoid Vektoren in DHDN
        //Querkrümmunsradius
        double N = aB / Math.sqrt(1.0 - e2B * Math.pow(Math.sin(B / 180 * Math.PI), 2));

        // Ergebnis Vektoren	
        xB = (N + h) * Math.cos(B / 180 * Math.PI) * Math.cos(L / 180 * Math.PI);
        yB = (N + h) * Math.cos(B / 180 * Math.PI) * Math.sin(L / 180 * Math.PI);
        zB = (N * (bB * bB) / (aB * aB) + h) * Math.sin(B / 180 * Math.PI);        
   
        //Umrechnung der Drehwinkel in Bogenmaß
        //http://de.wikipedia.org/wiki/Radiant_%28Einheit%29#Umrechnung_zwischen_Radiant_und_Grad
        double exRad = (ex * Math.PI / 180.0) / 3600.0;
        double eyRad = (ey * Math.PI / 180.0) / 3600.0;
        double ezRad = (ez * Math.PI / 180.0) / 3600.0;

        //Maßstabsumrechnung
        //http://de.wikipedia.org/wiki/Ma%C3%9Fstab_%28Verh%C3%A4ltnis%29
        double mEXP = 1 - m * Math.pow(10, -6);

        //Drehmatrix
        // 1         Ez    -Ez
        // -Ez       1      Ex 
        // Ey       -Ex     1

        //Rotierende Vektoren
        // = Drehmatrix * Vektoren in WGS84
        double RotVektor1 = 1.0 * xB + ezRad * yB + (-1.0 * eyRad * zB);
        double RotVektor2 = (-1.0 * ezRad) * xB + 1 * yB + exRad * zB;
        double RotVektor3 = (eyRad) * xB + (-1.0 * exRad) * yB + 1 * zB;    
        //Maßstab berücksichtigen
        double RotVectorM1 = RotVektor1 * mEXP;
        double RotVectorM2 = RotVektor2 * mEXP;
        double RotVectorM3 = RotVektor3 * mEXP;
        //Translation anbringen
        //dxT = Drehmatrix * dx * m
        double dxT = 1.0 * dx * mEXP + ezRad * dy * mEXP + (-1.0 * eyRad) * dz * mEXP;
        double dyT = (-1.0 * ezRad) * dx * mEXP + 1.0 * dy * mEXP + exRad * dz * mEXP;
        double dzT = (eyRad) * dx * mEXP + (-1.0 * exRad) * dy * mEXP + 1 * dz * mEXP;
        //Vektoren jetzt in WGS84
        xW = RotVectorM1 + dxT;
        yW = RotVectorM2 + dyT;
        zW = RotVectorM3 + dzT;
        
        double s = Math.sqrt(xW * xW + yW * yW);
        double T = Math.atan(zW * aW / (s * bW));
        double Bz = Math.atan((zW + e2W * (aW * aW) / bW * Math.pow(Math.sin(T), 3)) / (s - e2W * aW * Math.pow(Math.cos(T), 3)));

        double Lz = Math.atan(yW / xW);
        N = aW / Math.sqrt(1 - e2W * Math.pow(Math.sin(Bz), 2));

        Altitude = s / Math.cos(Bz);
        Latitude = Bz * 180 / Math.PI;
        Longitude = Lz * 180 / Math.PI;
  	
        
        log.debug("transformGK2WGS84():Longitude: "+Longitude);
        log.debug("transformGK2WGS84():Latitude: "+Latitude);
        log.debug("transformGK2WGS84():Bezugsmeridian: "+bMeridian);
        log.debug("transformGK2WGS84():Altitude: "+Altitude);
        
        returnArray.add(Longitude);
        returnArray.add(Latitude);
		return returnArray;
	}
	/***************************************************************************************************************/
	protected void finalize() throws Throwable {
	    try {
	        log = null;        //logger nullen, somit sollten auch File-Handles geschlossen werden
	        returnArray = null;
	    } finally {
	        super.finalize();
	    }
	}
}


