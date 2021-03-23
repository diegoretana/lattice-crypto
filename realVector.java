package ABB;

import java.io.Serializable;
import java.math.BigDecimal;

public class realVector implements Serializable
{
    private BigDecimal[] v;
    public int dim;
    
    public BigDecimal Norm2 = new BigDecimal("0");
    public double Norm;

    
    public realVector(BigDecimal[] x)
    {
        this.dim = x.length;
        this.v = x.clone();
        int i;
        
        for(i=0; i<dim; i++)
            Norm2 = Norm2.add( v[i].multiply(v[i]) );
        
        Norm = Math.sqrt( Norm2.doubleValue() );
    }
    
    public int getDim()
    {
        return dim;
    }

    public void print()
    {
        int i;
        int r = dim - 1;
        
        System.out.print("[");
        for(i=0; i<r; i++)
        {
            
            System.out.printf("\t%.4f", v[i].doubleValue() );
        }
        
        System.out.printf("\t%.4f ]", v[i].doubleValue() );
    }
    /*
    public LVector add(LVector u)
    {
        int i;

        long aux[] = new long[dim];

        for(i=0; i<dim; i++)
            aux[i] = (long)(v[i] + u.getAt(i));

        return new LVector(aux);
    }
    */
    public realVector add(realVector u)
    {
        int i;

        BigDecimal aux[] = new BigDecimal[dim];

        for(i=0; i<dim; i++)
            aux[i] = v[i].add( u.getAt(i) );

        return new realVector(aux);
    }
    
    public realVector sub(LVector u)
    {
        int i;

        BigDecimal aux[] = new BigDecimal[dim];

        for(i=0; i<dim; i++)
            aux[i] = v[i].subtract( new BigDecimal(u.getAt(i)) );

        return new realVector(aux);
    }
    
    
    public BigDecimal punto(LVector u)
    {
        int i;
        BigDecimal productopunto = new BigDecimal("0");
        
        if(u.getDim() == this.dim)
            for(i=0; i<dim; i++)
                productopunto = v[i].multiply( new BigDecimal(u.getAt(i)) );
        else
            System.err.println("no son de la misma dimension");
        
        return productopunto;
    }
    
    public BigDecimal punto(realVector u)
    {
        int i;
        BigDecimal productopunto = new BigDecimal("0");
        
        if(u.getDim() == this.dim)
            for(i=0; i<dim; i++)
                productopunto = productopunto.add( v[i].multiply( u.getAt(i) ) );
        else
            System.err.println("no son de la misma dimension");
        
        return productopunto;
    }
    
    /*
    public BigDecimal punto(Vector u)
    {
        int i;
        BigDecimal productopunto = new BigDecimal("0");
        
        if(u.getDim() == this.dim)
            for(i=0; i<dim; i++)
                productopunto = productopunto.add( v[i].multiply( u.getAt(i) ) );
        else
            System.err.println("no son de la misma dimension");
        
        return productopunto;
    }
    */
    
    
    public double getnorm()
    {
        BigDecimal norm2 = new BigDecimal("0");
        int i;
        
        for(i=0; i<dim; i++)
            norm2 = norm2.add( v[i].multiply( v[i] ) );
        
        Norm2 = norm2;
        Norm  = Math.sqrt( norm2.doubleValue() );
        
        return Norm;
    }

    
    public BigDecimal getAt(int i)
    {
        return v[i];
    }
    
     public realVector mult(BigDecimal k)
    {
        int i;
        
        BigDecimal aux[] = new BigDecimal[dim];

        for(i=0; i<dim; i++)
            aux[i] = v[i].multiply( k );
        
        return new realVector(aux);
    }
}