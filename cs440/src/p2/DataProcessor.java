package p2;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;
import java.util.Locale;

import javax.xml.bind.annotation.adapters.NormalizedStringAdapter;

public class DataProcessor {
    static String[] A1_cand = { "b", "a", "?" };
    static String[] A4_cand = { "u", "y", "l", "t", "?" };
    static String[] A5_cand = { "g", "p", "gg", "?" };
    static String[] A6_cand = { "c", "d", "cc", "i", "j", "k", "m", "r", "q",
            "w", "x", "e", "aa", "ff", "?" };
    static String[] A7_cand = { "v", "h", "bb", "j", "n", "z", "dd", "ff", "o", "?" };
    static String[] A9_cand = { "t", "f", "?" };
    static String[] A10_cand = { "t", "f", "?" };
    static String[] A12_cand = { "t", "f", "?" };
    static String[] A13_cand = { "g", "p", "s", "?" };

    public double[][] m_inputvs;
    public double[][] m_outputvs;

    private List<CreditData> m_datas;
    private List<LensData> m_lens;
    private List<BubilData> m_bubils;
    
    //Credit Data
    class CreditData {
    	public double[] m_inputs;
    	public double[] m_outputs;
    	
        public CreditData(double[] inputs, double[] outputs) {
            m_inputs = inputs;
            m_outputs = outputs;
        }
    }
    
    class BubilData{
    	
    	public double[] m_inputs;
    	public double[] m_outputs;
    	
        public BubilData(double[] inputs, double[] outputs) {
            m_inputs = inputs;
            m_outputs = outputs;
        }
    	
    }
    
    // Class for preprocessing Lens Data
    class LensData { 
    	public double[] m_inputs;
    	public double[] m_outputs;
    	public LensData(double[] inputs, double[] outputs) {
    		m_inputs = inputs;
    		m_outputs = outputs;
    	}
    }

    double cvtDouble(String [] candidates, String name) {
        for (int i = 0; i < candidates.length; ++i) {
            if (candidates[i].equals(name)) {
                return i;
            }
        }
        return candidates.length;
    }

    public DataProcessor(String aFileName, int id) throws FileNotFoundException {
    	
    	switch (id) {
    		// Test Data and Credit Data 
	    	default:
	    	case 0: 
		        m_datas = new ArrayList<CreditData>();
		        Scanner s0 = null;
		        
		        FileReader f0 = new FileReader(aFileName);
		        
		        try {
		            s0 = new Scanner(new BufferedReader(f0));
		            s0.useLocale(Locale.US);
		
		            while (s0.hasNextLine()) {
		                CreditData data = processCreditData(s0.nextLine());
		                m_datas.add(data);
		            }
		        } finally {
		            s0.close();
		        }
		
		        int i = 0;
		        m_inputvs = new double[m_datas.size()][];
		        m_outputvs = new double[m_datas.size()][];
		        for (CreditData data: m_datas) {
		            m_inputvs[i] = data.m_inputs;
		            m_outputvs[i] = data.m_outputs;
		            ++i;
		        }
		        m_inputvs = norm(m_inputvs);
		        break;
		        
		    // Lens Data  
	    	case 1: 
	    		m_lens = new ArrayList<LensData>();
	    		Scanner s1 = null;
	    		
	    		FileReader f1 = new FileReader(aFileName);
	    		
	    		try {
	    			s1 = new Scanner(new BufferedReader(f1));
	    			s1.useLocale(Locale.US);
	    			
	    			while (s1.hasNextLine()) {
	    				LensData dataLens = processLensData(s1.nextLine());
	    				m_lens.add(dataLens);
	    			}
	    		} finally {
	    			s1.close();
	    		}
	    		
	    		int j = 0;
	    		m_inputvs = new double[m_lens.size()][];
	    		m_outputvs = new double[m_lens.size()][];
	    		
	    		for (LensData data : m_lens) {
	    			m_inputvs[j] = data.m_inputs;
	    			m_outputvs[j] = data.m_outputs;
	    			++j;
	    		}
	    		m_inputvs = norm(m_inputvs);
	    		break;
	    		
	    		
	    		// Lens Data  
	    	case 2: 
	    		m_bubils = new ArrayList<BubilData>();
	    		Scanner s2 = null;
	    		
	    		FileReader f2 = new FileReader(aFileName);
	    		
	    		try {
	    			s2 = new Scanner(new BufferedReader(f2));
	    			s2.useLocale(Locale.US);
	    			
	    			while (s2.hasNextLine()) {
	    				BubilData dataBubils = processBubilData(s2.nextLine());
	    				m_bubils.add(dataBubils);
	    			}
	    		} finally {
	    			s2.close();
	    		}
	    		
	    		int k = 0;
	    		m_inputvs = new double[m_bubils.size()][];
	    		m_outputvs = new double[m_bubils.size()][];
	    		
	    		for (BubilData data : m_bubils) {
	    			m_inputvs[k] = data.m_inputs;
	    			m_outputvs[k] = data.m_outputs;
	    			++k;
	    		}
	    		m_inputvs = norm(m_inputvs);
	    		break;
    	}
    }

    private double nextDouble(Scanner s) {
        if (s.hasNextDouble()) {
            return s.nextDouble();
        } else {
            s.next();
            return 0.0;
        }
    }

    public CreditData processCreditData(String line) {
        Scanner scanner = new Scanner(line);
        scanner.useDelimiter(",");

        double[] inputs = new double[15];
        double[] outputs = new double[1];

        inputs[0] = cvtDouble(A1_cand, scanner.next());
        inputs[1] = nextDouble(scanner);
        inputs[2] = nextDouble(scanner);
        inputs[3] = cvtDouble(A4_cand, scanner.next());
        inputs[4] = cvtDouble(A5_cand, scanner.next());
        inputs[5] = cvtDouble(A6_cand, scanner.next());
        inputs[6] = cvtDouble(A7_cand, scanner.next());
        inputs[7] = nextDouble(scanner);
        inputs[8] = cvtDouble(A9_cand, scanner.next());
        inputs[9] = cvtDouble(A10_cand, scanner.next());
        inputs[10] = nextDouble(scanner);
        inputs[11] = cvtDouble(A12_cand, scanner.next());
        inputs[12] = cvtDouble(A13_cand, scanner.next());
        inputs[13] = nextDouble(scanner);
        inputs[14] = nextDouble(scanner);

        String output = scanner.next();
        if (output.equals("+")) {
            outputs[0] = 1.0;
        } else {
            outputs[0] = 0.0;
        }
        return new CreditData(inputs, outputs);
    }
    
    public LensData processLensData(String line) {
    	Scanner scanner = new Scanner(line);
    	scanner.useDelimiter(",");
    	
    	double[] inputs = new double[4];
    	double[] outputs = new double[1];
    	
    	for (int i = 0; i < 4; i++) {
			inputs[i] = scanner.nextDouble();
		}
    	
    	double output = scanner.nextDouble();
    	
    	if (output == 1) {
    		output = 0;
    	} else if (output == 2) {
    		output = 0.5;
    	} else {
    		output = 1;
    	}
    	
    	outputs[0] = output; 
    	return new LensData(inputs, outputs);
    }
    
    public BubilData processBubilData(String line) {
    	Scanner scanner = new Scanner(line);
    	scanner.useDelimiter(",");
    	
    	double[] inputs = new double[4];
    	double[] outputs = new double[1];
    	
    	for (int i = 0; i < 4; i++) {
			inputs[i] = scanner.nextDouble();
		}
    	
    	double output = scanner.nextDouble();
    	
    	if (output == 1) {
    		output = 0;
    	} else if (output == 2) {
    		output = 0.25;
    	} else if (output == 3) {
    		output = 0.5;
    	} else if (output == 4) {
    		output = 0.75;
    	} else {
    		output = 1;
    	}
    	
    	outputs[0] = output; 
    	return new BubilData(inputs, outputs);
    }
    
    // Normalizes double[][] inputvs
    public double[][] norm(double[][] inputvs) {
    	for (int i = 0; i < inputvs[0].length; i++) {
			double sum = 0;
			
			for (int j = 0; j < inputvs.length; j++) {
				sum += inputvs[j][i];
			}
			
			double n = 0;
			double avg = sum / inputvs.length; 
			
			for (int j = 0; j < inputvs.length; j++) {
				n += Math.pow(inputvs[j][i] - avg, 2);
			}
			
			double sigma = Math.sqrt(n / inputvs.length);
			
			for (int j = 0; j < inputvs.length; j++) {
				inputvs[j][i] = (inputvs[j][i] - avg) / sigma; 
			}
		}
    	return inputvs;
    }
}



