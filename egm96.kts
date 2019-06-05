import java.io.File
import java.io.InputStream
import kotlin.math.*

fun egm96(FLAT: Double, FLON: Double): Double {
val egm96 = File("egm96").reader().readLines()
val corrc = File("corrc").reader().readLines()

val NMAX = 360
val MM=((NMAX+1)*(NMAX+2))/2 +1
var COTHET: Double
var SITHET: Double
var GR: Double
var HACO: Double
var RE: Double
var X: Double
var Y: Double
var Z: Double
var FLATR: Double
var FLONR: Double
var RLAT: Double
var RLON: Double
var UNDU: Double
var AR: Double
var ARN: Double
var AC: Double
var A: Double
var SUM: Double
var SUMC: Double
var TEMP: Double
var TEMPC: Double
var T1: Double
var T2: Double
var N: Double

val RAD=57.29577951308232
val HT=0.0
val AE=6378137.0
val E2=.00669437999013
val GM=.3986004418E15
val GEQT=9.7803253359
val KK=.00193185265246

var P = DoubleArray(MM)
var HC = Array(MM) {0.0}
var HS = Array(MM) {0.0}
var CC = Array(MM) {0.0}
var CS = Array(MM) {0.0}
var RLEG = Array(NMAX+2) {0.0}
var RLNN = Array(NMAX+2) {0.0}
var SINML = Array(NMAX+2) {0.0}
var COSML = Array(NMAX+2) {0.0}
var DRTS = Array(1300) {0.0}
var DIRT = Array(1300) {0.0}
val lineRegex = "\\s+".toRegex()

for (line in corrc){
               val (n, m, t1, t2) = line.trim().split(lineRegex)
               val nInt = n.toInt()
               val mInt = m.toInt()
               val ig = (nInt * (nInt+ 1)) / 2 + mInt + 1
                    CC[ig] = t1.toDouble()
                    CS[ig] = t2.toDouble()
          }

for (line in egm96) {
                val (n, m, c1, c2) = line.trim().split(lineRegex)
                val nInt = n.toInt()
                val mInt = m.toInt()
                val ig = (nInt * (nInt + 1)) / 2 + mInt + 1
                HC[ig] = c1.toDouble()
                HS[ig] = c2.toDouble()
  }
      HC[4]=HC[4]+0.108262982131e-2/sqrt(5.0)
      HC[11]=HC[11]-.237091120053e-5/3.0
      HC[22]=HC[22]+0.608346498882e-8/sqrt(13.0)
      HC[37]=HC[37]-0.142681087920e-10/sqrt(17.0)
      HC[56]=HC[56]+0.121439275882e-13/sqrt(21.0)
      FLATR=FLAT/RAD
      FLONR=FLON/RAD
      T1= sin(FLATR).pow(2)
      N=AE/sqrt(1.0-E2*T1)
      T2=(N+HT)* cos(FLATR)
      X=T2*cos(FLONR)
      Y=T2*sin(FLONR)
      Z=(N*(1.0-E2)+HT)*sin(FLATR)
      RE=sqrt(X.pow(2)+Y.pow(2)+Z.pow(2))
      RLAT=atan(Z/sqrt(X.pow(2)+Y.pow(2)))
      GR=GEQT*(1.0+KK*T1)/ sqrt(1.0-E2*T1)
      RLAT=1.5707963267948966-RLAT

      for (n in 0..2*NMAX+1)
      {
      DRTS[n] = sqrt(n*1.0)
      DIRT[n] = 1.0/DRTS[n]
      }

      COTHET = cos(RLAT)
      SITHET = sin(RLAT)
      for (j in 0..NMAX) {
          RLNN[1] = 1.0
          RLNN[2] = SITHET * DRTS[3]

          for (n in 2..j) {
              RLNN[n + 1] = DRTS[2 * n + 1] * DIRT[2 * n] * SITHET * RLNN[n]
             }
              if (j + 2 <= NMAX + 1) {
                  RLEG[j + 1] = RLNN[j + 1]
              }
              if (j + 3 <= NMAX + 1) {
                  RLEG[j + 2] = DRTS[(j + 1) * 2 + 1] * COTHET * RLEG[j + 1]
              }
  for (m in j + 3..NMAX + 1) {
        RLEG[m]=DRTS[2*(m-1)+1]*DIRT[m+j-1]*DIRT[m-j-1]*(DRTS[2*(m-1)-1]*COTHET*RLEG[m-1]-DRTS[m+j-2]*DRTS[m-j-2]*DIRT[2*(m-1)-3]*RLEG[m-2])
         }
                  for (o in j + 1..NMAX + 1) {
                      P[(o - 1) * (o) / 2 + j + 1] = RLEG[o]
                  }

      }

      RLON=FLON/RAD
      SINML[1]=sin(RLON)
      COSML[1]= cos(RLON)
      SINML[2]=2.0*COSML[1]*SINML[1]
      COSML[2]=2.0*COSML[1]*COSML[1]-1.0

      for (p in 3..NMAX) {
          SINML[p] = 2.0 * COSML[1] * SINML[p - 1] - SINML[p - 2]
          COSML[p] = 2.0 * COSML[1] * COSML[p - 1] - COSML[p - 2]
      }

      AR=AE/RE
      ARN=AR
      AC=0.0
      A=0.0
      var K=3
      for (q in 2..NMAX) {
          ARN *= AR
          K += 1
          SUM = P[K] * HC[K]
          SUMC = P[K] * CC[K]
          for (r in 1..q) {
              K += 1
              TEMPC = CC[K] * COSML[r] + CS[K] * SINML[r]
              TEMP = HC[K] * COSML[r] + HS[K] * SINML[r]
              SUMC += P[K] * TEMPC
              SUM += P[K] * TEMP
          }

          AC += SUMC
          A += SUM * ARN
      }

      AC+=CC[1]+P[2]*CC[2]+P[3]*(CC[3]*COSML[1]+CS[3]*SINML[1])
      HACO=AC/100.0
      UNDU=A*GM/(GR*RE)
      UNDU+=HACO-0.530
      return UNDU
}

println(egm96(38.6281550,269.7791550))
