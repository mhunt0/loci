#ifndef BLOCK_HASH_H
#define BLOCK_HASH_H
#include <Tools/intervalSet.h>

namespace Loci {

  template <class T> class block_hash {

    struct block_info {
      unsigned long bits[16] ;
      T data[512] ;
      block_info() {
        for(int i=0;i<16;++i)
          bits[i] = 0 ;
      }
    };
    
    block_info **data[4096] ;
  
  public:
    block_hash() {
      for(int i=0;i<4096;++i)
        data[i] = 0 ;
    }
    void clear_hash() ;
      
    ~block_hash() {
      clear_hash() ;
    }

    void erase_set(intervalSet set) ;
    
    intervalSet domain() const ;

    T & access(int n) {
      const int a1 = (n) & 0x1ff ;
      const int a2 = (n>>9) &0x7ff ;
      const int a3 = (n>>20) &0xfff ;
      if(data[a3]!=0 && data[a3][a2]!=0) {
        data[a3][a2]->bits[(a1>>5)&0xf] |= (1<<(a1&0x1f)) ;
        return data[a3][a2]->data[a1] ;
      }
      if(data[a3] == 0) {
        data[a3] = new block_info *[2048] ;
        for(int i=0;i<2048;++i)
          data[a3][i] = 0 ;
      }
      if(data[a3][a2] == 0) 
        data[a3][a2] = new block_info ;
      data[a3][a2]->bits[(a1>>5)&0xf] |= (1<<(a1&0x1f)) ;
      return data[a3][a2]->data[a1] ;
    }

    void copy_hash(const block_hash &cp) ;
    block_hash(const block_hash &cp) {
      for(int i=0;i<4096;++i)
        data[i] = 0 ;
      copy_hash(cp) ;
    }
    block_hash &operator=(block_hash &cp) {
      copy_hash(cp) ;
    }

    T & elem(int n) {
      const int a1 = (n) & 0x1ff ;
      const int a2 = (n>>9) &0x7ff ;
      const int a3 = (n>>20) &0xfff ;

      return data[a3][a2]->data[a1] ;
    }

    const T & elem(int n) const {
      const int a1 = (n) & 0x1ff ;
      const int a2 = (n>>9) &0x7ff ;
      const int a3 = (n>>20) &0xfff ;

      return data[a3][a2]->data[a1] ;
    }
    const T &operator[](int indx) const {return elem(indx) ;}
    T &operator[](int indx) {return access(indx) ;}

    const T &operator()(int indx) const {return elem(indx) ;}
    T &operator()(int indx) {return elem(indx);}
  } ;

  template <class T> void block_hash<T>::clear_hash() {
    for(int i=0;i<4096;++i) {
      if(data[i]!=0) {
        for(int j=0;j<2048;++j) {
          if(data[i][j]!=0) {
            delete data[i][j] ;
          }
        }
        delete[] data[i] ;
      }
      data[i] = 0 ;
    }
  }

  template <class T> void block_hash<T>::copy_hash(const block_hash<T> &cp) {
    clear_hash() ;
    for(int i=0;i<4096;++i) {
      if(cp->data[i] != 0) {
        data[i] = new block_info *[2048] ;
        for(int j=0;j<2048;++j) {
          if(cp->data[i][j] != 0)
            data[i][j] = new block_info(cp->data[i][j]) ;
        }
      }
    }
  }
  template <class T> void block_hash<T>::erase_set(intervalSet set) {
    for(intervalSet::const_iterator ei = set.begin();
        ei != set.end() ;
        ++ei) {
      const int a1 = (*ei) & 0x1ff ;
      const int a2 = (*ei>>9) &0x7ff ;
      const int a3 = (*ei>>20) &0xfff ;
      data[a3][a2]->bits[(a1>>5)&0xf] &= ~(1<<(a1&0x1f)) ;
    }
    // Clean up any fully erased blocks
    for(int i=0;i<4096;++i) {
      if(data[i]!=0) {
        for(int j=0;j<2048;++j) {
          if(data[i][j]!=0) {
            bool all_zero = true ;
            for(int k=0;k<16;++k)
              if(data[i][j]->bits[k] != 0)
                all_zero = false ;
            if(all_zero) {
              delete data[i][j] ;
              data[i][j] = 0 ;
            }
          }
        }
      }
    }
      
  }    
  
  template <class T> intervalSet block_hash<T>::domain() const  {
    intervalSet val ;

    for(int i=0;i<4096;++i) {
      if(data[i]!=0) {
        int a3 = i << 20 ;
        for(int j=0;j<2048;++j) {
          if(data[i][j]!=0) {
            int a2 = j << 9 ;
            for(int k=0;k<16;++k) {
              const unsigned long bit = data[i][j]->bits[k] ;
              if(bit == 0) 
                continue ;
              int a1 = k<<5 ;
              if(bit == 0xffffffff) {
                int low = a1|a2|a3 ;
                val += interval(low,low+32-1) ;
              } else {
                for(int l=0;l<32;++l) {
                  if((bit & (1<<l)) != 0)
                    val += (a1|a2|a3|l) ;
                }
              }
            }
          }
        }
      }
    }
    return val ;
  }
                 
                

}
#endif
