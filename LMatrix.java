package ABB;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.StringTokenizer;
import java.math.BigDecimal;
import java.math.RoundingMode;

public class LMatrix implements Serializable
{
    public  long[][] M;
    public  int rows;
    public  int columns;
    
    public int fill = 4;
    public realVector[] gso;
    public LVector[] basis;
    
    private static int scale = SystemParameters.scale;
    
    public LMatrix()
    {
        
    }
    
    public LMatrix(int i, int j)
    {
        M = new long[i][j];
    }
    
    public LMatrix(long[][] mat)
    {
        M = mat;
        rows = mat.length;
        columns = mat[0].length;
    }
    
    private long Mod(long a, long q)
    {
        long mod;
        mod = a % q;
        mod = mod >= 0 ? mod : q + mod;
        
        return mod;
    }
    
    public void setIdentity()
    {
        int i, j;
        BigDecimal[] gsovec;
        this.gso = new realVector[rows];
        
        for(i=0; i<rows; i++)
        {
            M[i][i] = 1;
            gsovec = new BigDecimal[rows];
            for(j=0; j<rows; j++)
                gsovec[j] = new BigDecimal("0");
            gsovec[i] = new BigDecimal("1");
            this.gso[i] = new realVector( gsovec );
        }
        
        toBasis();
    }
    
    public LMatrix addition(LMatrix B)
    {
        long[][] C = new long[rows][columns];
        
        int i, j;
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                C[i][j] = M[i][j] + B.M[i][j];
        
        return new LMatrix(C);
    }
    
    public void mult_const(long k)
    {        
        for(int i=0; i<rows; i++)
            for(int j=0; j<columns; j++)
                M[i][j] = k*M[i][j];
    }
    
    public LMatrix mult(LMatrix B)
    {
        int r = rows;
        int c = B.columns;
        
        long C[][] = new long[r][c];
        
        int i, j, k;
        for(i=0; i<r; i++)
            for(j=0; j<c; j++)
                for(k=0; k<columns; k++)
                    C[i][j] += M[i][k] * B.M[k][j];
        
        return new LMatrix(C);
    }
    
    public LVector ModMult(LVector v)
    {
        long[] u = new long[rows];
        long q = SystemParameters.q;
        
        int i, j;
        for(i=0; i<rows; i++)
        {
            for(j=0; j<columns; j++)
                u[i] += M[i][j] * v.getAt(j);
            
            u[i] = Mod( u[i], q);
        }
        
        return new LVector(u);
    }
    
    public LVector Mult(LVector v)
    {
        long[] u = new long[rows];
        
        int i, j;
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                u[i] += M[i][j] * v.getAt(j);
        
        return new LVector(u);
    }
    
    
    public double gsoNorm()
    {
        double maxnorm = 0, auxnorm;
        int i;
        
        if(gso != null)
            for(i=0; i<gso.length; i++)
            {
                auxnorm = gso[i].getnorm();
                if(auxnorm > maxnorm)
                    maxnorm = auxnorm;
            }
        
        return maxnorm;
    }

    
    public LMatrix ModMult(LMatrix B, long modulo)
    {
        int r = rows;
        int c = B.columns;
        
        long C[][] = new long[r][c];
        
        int i, j, k;
        for(i=0; i<r; i++)
            for(j=0; j<c; j++)
            {
                for(k=0; k<columns; k++)
                    C[i][j] += M[i][k] * B.M[k][j];
                
                C[i][j] = Mod( C[i][j], modulo );
            }
        
        return new LMatrix(C);
    }
    
    public LMatrix augment(LMatrix B)
    {
        long[][] C = new long[rows][columns + B.columns];
        
        int i, j;
        
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                C[i][j] = M[i][j];
        
        for(i=0; i<rows; i++)
            for(j=0; j<B.columns; j++)
                C[i][columns + j] = B.M[i][j];
        
        return new LMatrix(C);
    }
    
    public LMatrix augmentDown(LMatrix B)
    {
        long[][] C = new long[rows + B.rows][columns];
        
        int i, j;
        
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                C[i][j] = M[i][j];
        
        for(i=0; i<B.rows; i++)
            for(j=0; j<B.columns; j++)
                C[rows + i][j] = B.M[i][j];
        
        return new LMatrix(C);
    }    
    
    public LMatrix transpose()
    {
        long[][] T = new long[columns][rows];
        
        int i, j;
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                T[j][i] = M[i][j];
        
        return new LMatrix(T);
    }
    
    public void toBasis()
    {
        this.basis = new LVector[columns];
        //double maxnorm = 0, aux;
        
        long[] v;
        
        int i, j;
        for(j=0; j<columns; j++)
        {
            v = new long[rows];
            for(i=0; i<rows; i++)
                v[i] = M[i][j];
            
            basis[j] = new LVector(v);
            //aux = basis[j].norm();
            //if( aux > maxnorm )
            //    maxnorm = aux;
        }
        
        //System.out.println(" || Ta || = " + String.valueOf( maxnorm ) );
    }
    
    public void Compute_Gram_Shmidt()
    {
        int i, j;
        BigDecimal muij;
        
        toBasis();
        
        gso = new realVector[columns];
        
        gso[0] = basis[0].torealVector();
        
        realVector aux;
        for(i=1; i<columns; i++)
        {
            muij = gso[0].punto( basis[i] ).divide( gso[0].Norm2, scale, RoundingMode.HALF_UP );
            aux = gso[0].mult( muij );
            for(j=1; j<i; j++)
            {
                muij = gso[j].punto( basis[i] ).divide( gso[j].Norm2, scale, RoundingMode.HALF_UP );
                aux = aux.add( gso[j].mult(muij) );
            }
            
            gso[i] = basis[i].sub(aux);
            System.out.println(" gso . . . " + i + " / " + String.valueOf(gso.length) );
        }
    }
    
    public String getString()
    {
        String str = "{";
        
        int i, j;
        
        for(i=0; i<rows; i++)
        {
            str = str + "{";
            for(j=0; j<columns; j++)
                str = str + String.valueOf(M[i][j]) + ", ";
            
            str = str.substring(0, str.length()-2);
            str = str + "}, ";
        }
        str = str.substring(0, str.length()-2);
        str = str + "}";
        
        return str;
    }
    
    public void WriteLMatrix(String path)
    {
        int i, j, last=rows-1;
        
        File file = new File(path);
        FileWriter writer = null;
        
        String str = "  " + String.valueOf(rows) + "  " + String.valueOf(columns) + " \n{";
        
        try
        {
            file.createNewFile();
            writer = new FileWriter(file);
            writer.write(str);
            writer.flush();
        }catch(IOException ioe){System.err.println(" Error : " + ioe.getMessage());}
        
        for(i=0; i<rows; i++)
        {
            str = "{";
            for(j=0; j<columns; j++)
                str = str + String.valueOf(M[i][j]) + ", ";
            
            if( i == last )
            {
                str = str.substring(0, str.length()-2);
                str = str + "}}";
            }
            else
            {
                str = str.substring(0, str.length()-2);
                str = str + "}, ";
            }
            
            try{
                writer.write(str);
                writer.flush();
            }catch(IOException ioe){System.err.println(" Error : " + ioe.getMessage());}
        }
        
        try{
            writer.close();
        }catch(IOException ioe){System.err.println(" Error : " + ioe.getMessage());}
    }
    
    public void print()
    {
        int i, j;
        long aux;
        String cad;

        for(i=0; i<rows; i++)
        {
            for(j=0; j<columns; j++)
            {
                aux = M[i][j];
                if(aux == 0)
                    cad = "_";
                else
                    cad = String.valueOf(aux);
                while(cad.length() < fill)
                    cad = " " + cad;
                System.out.print(cad);
            }
            System.out.println("");
        }
    }
    
    public void printGSO()
    {
        int i;
        
        for(i=0; i<gso.length; i++)
        {
            System.out.printf(" %d : %.4f : ", i, gso[i].getnorm() );
            gso[i].print();
            System.out.println("");
        }
    }
    
    public void reduceModule(int q)
    {
        int i, j;
        
        for(i=0; i<rows; i++)
            for(j=0; j<columns; j++)
                M[i][j] = Mod( M[i][j], q );
    }
    
    public static long[][] readMatrix(String path)
    {
        File fichero = new File(path);
        StringTokenizer token;
        String line, s;
        int filas=0, columnas=0, fila=0, col=0, max=0;
        long num;
        long[][] matriz = null;     
        
        try
        {
            FileReader fr = new FileReader(fichero);
            BufferedReader br = new BufferedReader(fr);
            
            line = br.readLine();
            token = new StringTokenizer(line);
            filas = Integer.parseInt(token.nextToken());
            columnas = Integer.parseInt(token.nextToken());
            matriz = new long[filas][columnas];
            
            while( (line = br.readLine()) != null)
            {
                token = new StringTokenizer(line," ,{}");
                
                while( token.hasMoreElements() )
                {
                    s = token.nextToken();
                    if( s.length() > max)
                        max = s.length();
                    
                    num = Long.parseLong( s );
                    matriz[fila][col] = num;
                    
                    col = (++col)%columnas;
                    if(col == 0)
                        fila++;
                }
            }            
        }catch(IOException ioe)
        {
            System.err.println( ioe.getMessage() );
        }
        
        return matriz;
    }
    
    public static void main(String[] args)
    {
        long[][] m1 = 
        {   {183, 1, 126, 26, 96, 115}, 
            {49, 241, 130, 130, 154, 111},
            {154, 20, 243, 257, 11, 4}, 
            {92, 96, 91, 131, 89, 136}, 
            {232, 100, 14, 119, 63, 23}, 
            {95, 22, 100, 82, 108, 15}};
        
        long[][] iTa = readMatrix("C:\\Users\\Diego\\Dropbox\\alpha\\Lattice\\IBE-Workspace\\01_3_504_3971791\\S_01_504_504_3971791.tdr");
        
        LMatrix mat1 = new LMatrix(iTa);
        LMatrix aux = mat1.transpose();
        
        mat1.Compute_Gram_Shmidt();
        
        aux.print(); System.out.println("");
        mat1.printGSO();
        
    }

}
