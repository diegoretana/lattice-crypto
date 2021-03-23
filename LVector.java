package ABB;

import java.math.BigDecimal;
import java.io.Serializable;

public class LVector implements Serializable
{
    private long v[];
    private int dim;

    
    public LVector(int dim)
    {
        this.dim = dim;
        v = new long[dim];
    }
    
    public LVector(long[] x)
    {
        this.dim = x.length;
        this.v = x.clone();
    }
    
    private long Mod(long a, long q)
    {
        long mod;
        mod = a % q;
        mod = mod >= 0 ? mod : q + mod;
        
        return mod;
    }
    
    public int getDim()
    {
        return dim;
    }

    public void print()
    {
        int tabulador = 7;
        int i, p;
        int r = dim - 1;
        String s;
        
        System.out.print("{");
        for(i=0; i<r; i++)
        {
            p = (int)v[i];
            s = String.valueOf(p) + ",";
            while(s.length() < tabulador)
                s = " " + s;
            System.out.print(s);
        }
        p = (int)v[i];
        s = String.valueOf(p);
        while(s.length() < tabulador)
            s = " " + s;
        System.out.printf(s + "}");
    }
    
    public LVector sub(LVector u)
    {
        long a, q = SystemParameters.q;
        int i;
        LVector vaux = null;
        
        if(u.dim == this.dim)
        {
            long aux[] = new long[dim];
            
            for(i=0; i<dim; i++)
            {
                a = v[i] - u.getAt(i);
                aux[i] = Mod(a, q);
            }
            
            vaux = new LVector(aux);
        }
        else
        {
            System.err.println("no son de la misma dimension");
        }
        
        return vaux;
    }
    
    public realVector sub(realVector u)
    {
        int i;
        
        BigDecimal aux[] = new BigDecimal[dim];
        BigDecimal ui, _1 = new BigDecimal("-1");

        for(i=0; i<dim; i++)
        {
            ui = u.getAt(i).multiply( _1 );
            aux[i] = ui.add( new BigDecimal( v[i] ) );
        }

        return new realVector(aux);
    }
    
    public LVector add(LVector u)
    {
        int i;
        LVector vaux = null;
        
        if(u.dim == this.dim)
        {
            long aux[] = new long[dim];
            
            for(i=0; i<dim; i++)
                aux[i] = v[i] + u.getAt(i);
            
            vaux = new LVector(aux);
        }
        else
        {
            System.err.println("no son de la misma dimension");
        }
        
        return vaux;
    }
    
    public LVector mult(int k)
    {
        int i;
        LVector vaux;
        long aux[] = new long[dim];

        for(i=0; i<dim; i++)
            aux[i] = k*v[i];

        vaux = new LVector(aux);
        
        return vaux;
    }
    
    public LVector mult(long k)
    {
        int i;
        LVector vaux;
        long aux[] = new long[dim];

        for(i=0; i<dim; i++)
            aux[i] = k*v[i];

        vaux = new LVector(aux);
        
        return vaux;
    }
    
    public LVector augment(LVector x)
    {
        int len = this.dim + x.getDim();
        
        long[] vx = new long[len];
        
        int i, k;
        for(i=0, k=0; i<this.dim; i++, k++)
            vx[k] = v[i];
        for(i=0; i<x.getDim(); i++, k++)
            vx[k] = x.getAt(i);
        
        return new LVector(vx);
    }
    
    public long norm2()
    {
        long norm2 = 0;
        int i;
        
        for(i=0; i<dim; i++)
            norm2 += v[i]*v[i];
        
        return norm2;
    }
    
    public double norm()
    {
        return Math.sqrt(  norm2() );
    }
    
    public long punto(LVector u)
    {
        int i;
        long productopunto = 0;
        
        if(u.dim == this.dim)
            for(i=0; i<dim; i++)
                productopunto += v[i]*u.getAt(i);
        else
            System.err.println("no son de la misma dimension");
        
        return productopunto;
    }
    
    public long modpunto(LVector u)
    {
        int i;
        long productopunto = 0;
        long q = SystemParameters.q;
        
        if(u.dim == this.dim)
            for(i=0; i<dim; i++)
            {
                productopunto += v[i]*u.getAt(i);
                productopunto = Mod(productopunto,q);
            }
        else
            System.err.println("no son de la misma dimension");
        
        return productopunto;
    }
    
    public long getAt(int i)
    {
        return v[i];
    }
    
    public void reduceModule(long q)
    {
        long aux;
        int i;
        
        for(i=0; i<dim; i++)
            v[i] = Mod(v[i],q);
    }
    
    public realVector torealVector()
    {
        BigDecimal[] d = new BigDecimal[dim];
        
        for(int i=0; i<dim; i++)
            d[i] = new BigDecimal( v[i] );
        
        return new realVector(d);
    }
    
    public boolean iszero()
    {
        int i;
        double p=1;
        
        for(i=0; i<dim; i++)
            p *= v[i];
        
        if(p == 0)
            return true;
        else
            return false;
    }
    
    public String getString()
    {
        String strvec = "{";
        int i;
        for(i=0; i<v.length; i++)
            strvec += String.valueOf(v[i]) + ",";
        strvec = strvec.substring(0, strvec.length()-1);
        strvec += "}";
        
        return strvec;
    }
    
}
