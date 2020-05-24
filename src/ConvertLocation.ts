import GeoPoint from './GeoPoint';

export default class ConvertLocation {
  public GEO = 0;
  public KATEC = 1;
  public TM = 2;
  public GRS80 = 3;
  public UTMK = 4;

  private mInd: number[] = [];
  // tslint:disable-next-line: variable-name
  private m_Es: number[] = [];
  // tslint:disable-next-line: variable-name
  private m_Esp: number[] = [];
  private srcM: number[] = [];
  private dstM: number[] = [];

  private EPSLN = 0.0000000001;
  private mArMajor: number[] = [];
  private mArMinor: number[] = [];

  private mArScaleFactor: number[] = [];
  private mArLonCenter: number[] = [];
  private mArLatCenter: number[] = [];
  private mArFalseNorthing: number[] = [];
  private mArFalseEasting: number[] = [];

  private datumParams: number[] = [];

  constructor() {
    // GEO Initial
    this.mArScaleFactor[this.GEO] = 1;
    this.mArLonCenter[this.GEO] = 0.0;
    this.mArLatCenter[this.GEO] = 0.0;
    this.mArFalseNorthing[this.GEO] = 0.0;
    this.mArFalseEasting[this.GEO] = 0.0;
    this.mArMajor[this.GEO] = 6378137.0;
    this.mArMinor[this.GEO] = 6356752.3142;

    // KATEC Initial
    this.mArScaleFactor[this.KATEC] = 0.9999;
    this.mArLonCenter[this.KATEC] = 2.23402144255274; // 128
    this.mArLatCenter[this.KATEC] = 0.663225115757845;
    this.mArFalseNorthing[this.KATEC] = 600000.0;
    this.mArFalseEasting[this.KATEC] = 400000.0;
    this.mArMajor[this.KATEC] = 6377397.155;
    this.mArMinor[this.KATEC] = 6356078.9633422494;

    // TM Initial
    this.mArScaleFactor[this.TM] = 1.0;
    this.mArLonCenter[this.TM] = 2.21661859489671;
    this.mArLatCenter[this.TM] = 0.663225115757845;
    this.mArFalseNorthing[this.TM] = 500000.0;
    this.mArFalseEasting[this.TM] = 200000.0;
    this.mArMajor[this.TM] = 6377397.155;
    this.mArMinor[this.TM] = 6356078.9633422494;

    // GRS80 Initial
    this.mArScaleFactor[this.GRS80] = 1.0;
    this.mArLonCenter[this.GRS80] = 2.2165681500328;
    this.mArLatCenter[this.GRS80] = 0.663225115757845;
    this.mArFalseNorthing[this.GRS80] = 500000.0;
    this.mArFalseEasting[this.GRS80] = 200000.0;
    this.mArMajor[this.GRS80] = 6378137;
    this.mArMinor[this.GRS80] = 6356752.3142;

    // UTMK Initial
    this.mArScaleFactor[this.UTMK] = 0.9996;
    this.mArLonCenter[this.UTMK] = 2.22529479629277;
    this.mArLatCenter[this.UTMK] = 0.663225115757845;
    this.mArFalseNorthing[this.UTMK] = 2000000.0;
    this.mArFalseEasting[this.UTMK] = 1000000.0;
    this.mArMajor[this.UTMK] = 6378137;
    this.mArMinor[this.UTMK] = 6356752.3141403558;

    this.datumParams[0] = -146.43;
    this.datumParams[1] = 507.89;
    this.datumParams[2] = 681.46;

    this.initialEsEspInd(this.GEO);
    this.initialEsEspInd(this.KATEC);
    this.initialEsEspInd(this.TM);
    this.initialEsEspInd(this.UTMK);
    this.initialEsEspInd(this.GRS80);

    this.initialSrcDst(this.GEO);
    this.initialSrcDst(this.KATEC);
    this.initialSrcDst(this.TM);
    this.initialSrcDst(this.UTMK);
    this.initialSrcDst(this.GRS80);
  }

  private initialEsEspInd(type: number) {
    const tmp = this.mArMinor[type] / this.mArMajor[type];

    this.m_Es[type] = 1.0 - tmp * tmp;
    this.m_Esp[type] = this.m_Es[type] / (1.0 - this.m_Es[type]);

    if (this.m_Es[type] < 0.00001) {
      this.mInd[type] = 1.0;
    } else {
      this.mInd[type] = 0.0;
    }
  }

  private initialSrcDst(type: number) {
    const es = this.m_Es[type];
    const data =
      this.mArMajor[type] *
      this.mlfn(this.e0fn(es), this.e1fn(es), this.e2fn(es), this.e3fn(es), this.mArLatCenter[type]);
    this.srcM[type] = data;
    this.dstM[type] = data;
  }

  private D2R(degree: number) {
    return (degree * Math.PI) / 180.0;
  }

  private R2D(radian: number) {
    return (radian * 180.0) / Math.PI;
  }

  private mlfn(e0: number, e1: number, e2: number, e3: number, phi: number) {
    return e0 * phi - e1 * Math.sin(2.0 * phi) + e2 * Math.sin(4.0 * phi) - e3 * Math.sin(6.0 * phi);
  }

  private e0fn(x: number) {
    return 1.0 - 0.25 * x * (1.0 + (x / 16.0) * (3.0 + 1.25 * x));
  }

  private e1fn(x: number) {
    return 0.375 * x * (1.0 + 0.25 * x * (1.0 + 0.46875 * x));
  }

  private e2fn(x: number) {
    return 0.05859375 * x * x * (1.0 + 0.75 * x);
  }

  private e3fn(x: number) {
    return x * x * x * (35.0 / 3072.0);
  }

  private asinz(value: number) {
    if (Math.abs(value) > 1.0) value = value > 0 ? 1 : -1;
    return Math.asin(value);
  }

  public convert(srctype: number, dsttype: number, inPt: GeoPoint) {
    const tmpPt = new GeoPoint();
    const outPt = new GeoPoint();

    if (srctype === this.GEO) {
      tmpPt.x = this.D2R(inPt.x);
      tmpPt.y = this.D2R(inPt.y);
    } else {
      this.tm2geo(srctype, inPt, tmpPt);
    }

    if (dsttype === this.GEO) {
      outPt.x = this.R2D(tmpPt.x);
      outPt.y = this.R2D(tmpPt.y);
    } else {
      this.geo2tm(dsttype, tmpPt, outPt);
    }

    return outPt;
  }

  public geo2tm(dsttype: number, inPt: GeoPoint, outPt: GeoPoint) {
    let x: number;
    let y: number;

    this.transform(this.GEO, dsttype, inPt);

    const deltaLon = inPt.x - this.mArLonCenter[dsttype];
    const sinPhi = Math.sin(inPt.y);
    const cosPhi = Math.cos(inPt.y);

    if (this.mInd[dsttype] !== 0) {
      const b = cosPhi * Math.sin(deltaLon);

      if (Math.abs(Math.abs(b) - 1.0) < this.EPSLN) {
        // 무한대 에러
        // tslint:disable-next-line: no-console
        console.error('무한대 에러');
      }
    } else {
      const b = 0;
      x = 0.5 * this.mArMajor[dsttype] * this.mArScaleFactor[dsttype] * Math.log((1.0 + b) / (1.0 - b));

      // tslint:disable-next-line: no-shadowed-variable
      let con = Math.acos((cosPhi * Math.cos(deltaLon)) / Math.sqrt(1.0 - b * b));

      if (inPt.y < 0) {
        con = con * -1;
        y = this.mArMajor[dsttype] * this.mArScaleFactor[dsttype] * (con - this.mArLatCenter[dsttype]);
      }
    }

    const al = cosPhi * deltaLon;
    const als = al * al;
    const c = this.m_Es[dsttype] * cosPhi * cosPhi;
    const tq = Math.tan(inPt.y);
    const t = tq * tq;
    const con = 1.0 - this.m_Es[dsttype] * sinPhi * sinPhi;
    const n = this.mArMajor[dsttype] / Math.sqrt(con);
    const ml =
      this.mArMajor[dsttype] *
      this.mlfn(
        this.e0fn(this.m_Es[dsttype]),
        this.e1fn(this.m_Es[dsttype]),
        this.e2fn(this.m_Es[dsttype]),
        this.e3fn(this.m_Es[dsttype]),
        inPt.y,
      );

    outPt.x =
      this.mArScaleFactor[dsttype] *
        n *
        al *
        (1.0 +
          (als / 6.0) *
            (1.0 - t + c + (als / 20.0) * (5.0 - 18.0 * t + t * t + 72.0 * c - 58.0 * this.m_Esp[dsttype]))) +
      this.mArFalseEasting[dsttype];

    outPt.y =
      this.mArScaleFactor[dsttype] *
        (ml -
          this.dstM[dsttype] +
          n *
            tq *
            (als *
              (0.5 +
                (als / 24.0) *
                  (5.0 -
                    t +
                    9.0 * c +
                    4.0 * c * c +
                    (als / 30.0) * (61.0 - 58.0 * t + t * t + 600.0 * c - 330.0 * this.m_Esp[dsttype]))))) +
      this.mArFalseNorthing[dsttype];
  }

  public tm2geo(srctype: number, inPt: GeoPoint, outPt: GeoPoint) {
    const tmpPt = new GeoPoint(inPt.x, inPt.y);
    const maxIter = 6;

    if (this.mInd[srctype] !== 0) {
      const f = Math.exp(inPt.x / (this.mArMajor[srctype] * this.mArScaleFactor[srctype]));
      const g = 0.5 * (f - 1.0 / f);
      const temp = this.mArLatCenter[srctype] + tmpPt.y / (this.mArMajor[srctype] * this.mArScaleFactor[srctype]);
      const h = Math.cos(temp);
      // tslint:disable-next-line: no-shadowed-variable
      const con = Math.sqrt((1.0 - h * h) / (1.0 + g * g));
      outPt.y = this.asinz(con);

      if (temp < 0) outPt.y *= -1;

      if (g === 0 && h === 0) {
        outPt.x = this.mArLonCenter[srctype];
      } else {
        outPt.x = Math.atan(g / h) + this.mArLonCenter[srctype];
      }
    }

    tmpPt.x -= this.mArFalseEasting[srctype];
    tmpPt.y -= this.mArFalseNorthing[srctype];

    const con = (this.srcM[srctype] + tmpPt.y / this.mArScaleFactor[srctype]) / this.mArMajor[srctype];

    let phi = con;
    let i = 0;

    while (true) {
      const delatPhi =
        (con +
          this.e1fn(this.m_Es[srctype]) * Math.sin(2.0 * phi) -
          this.e2fn(this.m_Es[srctype]) * Math.sin(4.0 * phi) +
          this.e3fn(this.m_Es[srctype]) * Math.sin(6.0 * phi)) /
          this.e0fn(this.m_Es[srctype]) -
        phi;

      phi = phi + delatPhi;

      if (Math.abs(delatPhi) <= this.EPSLN) break;

      if (i >= maxIter) {
        break;
      }

      i++;
    }

    if (Math.abs(phi) < Math.PI / 2) {
      const sinPhi = Math.sin(phi);
      const cosPhi = Math.cos(phi);
      const tanPhi = Math.tan(phi);
      const c = this.m_Esp[srctype] * cosPhi * cosPhi;
      const cs = c * c;
      const t = tanPhi * tanPhi;
      const ts = t * t;
      const cont = 1.0 - this.m_Es[srctype] * sinPhi * sinPhi;
      const n = this.mArMajor[srctype] / Math.sqrt(cont);
      const r = (n * (1.0 - this.m_Es[srctype])) / cont;
      const d = tmpPt.x / (n * this.mArScaleFactor[srctype]);
      const ds = d * d;

      outPt.y =
        phi -
        ((n * tanPhi * ds) / r) *
          (0.5 -
            (ds / 24.0) *
              (5.0 +
                3.0 * t +
                10.0 * c -
                4.0 * cs -
                9.0 * this.m_Esp[srctype] -
                (ds / 30.0) * (61.0 + 90.0 * t + 298.0 * c + 45.0 * ts - 252.0 * this.m_Esp[srctype] - 3.0 * cs)));
      outPt.x =
        this.mArLonCenter[srctype] +
        (d *
          (1.0 -
            (ds / 6.0) *
              (1.0 +
                2.0 * t +
                c -
                (ds / 20.0) * (5.0 - 2.0 * c + 28.0 * t - 3.0 * cs + 8.0 * this.m_Esp[srctype] + 24.0 * ts)))) /
          cosPhi;
    } else {
      outPt.y = Math.PI * 0.5 * Math.sin(tmpPt.y);
      outPt.x = this.mArLonCenter[srctype];
    }

    this.transform(srctype, this.GEO, outPt);
  }

  public getDistancebyGeo(pt1: GeoPoint, pt2: GeoPoint) {
    const lat1 = this.D2R(pt1.y);
    const lon1 = this.D2R(pt1.x);
    const lat2 = this.D2R(pt2.y);
    const lon2 = this.D2R(pt2.x);

    const longitude = lon2 - lon1;
    const latitude = lat2 - lat1;

    const a =
      Math.pow(Math.sin(latitude / 2.0), 2) + Math.cos(lat1) * Math.cos(lat2) * Math.pow(Math.sin(longitude / 2.0), 2);

    return 6376.5 * 2.0 * Math.atan2(Math.sqrt(a), Math.sqrt(1.0 - a));
  }

  public getDistancebyKatec(pt1: GeoPoint, pt2: GeoPoint) {
    pt1 = this.convert(this.KATEC, this.GEO, pt1);
    pt2 = this.convert(this.KATEC, this.GEO, pt2);

    return this.getDistancebyGeo(pt1, pt2);
  }

  public getDistancebyTm(pt1: GeoPoint, pt2: GeoPoint) {
    pt1 = this.convert(this.TM, this.GEO, pt1);
    pt2 = this.convert(this.TM, this.GEO, pt2);

    return this.getDistancebyGeo(pt1, pt2);
  }

  public getDistancebyUTMK(pt1: GeoPoint, pt2: GeoPoint) {
    pt1 = this.convert(this.UTMK, this.GEO, pt1);
    pt2 = this.convert(this.UTMK, this.GEO, pt2);

    return this.getDistancebyGeo(pt1, pt2);
  }

  public getDistancebyGrs80(pt1: GeoPoint, pt2: GeoPoint) {
    pt1 = this.convert(this.GRS80, this.GEO, pt1);
    pt2 = this.convert(this.GRS80, this.GEO, pt2);

    return this.getDistancebyGeo(pt1, pt2);
  }

  private getTimebySec(distance: number) {
    return Math.round((3600 * distance) / 4);
  }

  public getTimebyMin(distance: number) {
    return Math.ceil(this.getTimebySec(distance) / 60);
  }

  /*
	Author:       Richard Greenwood rich@greenwoodmap.com
	License:      LGPL as per: http://www.gnu.org/copyleft/lesser.html
	*/

  /**
   * convert between geodetic coordinates (longitude, latitude, height)
   * and gecentric coordinates (X, Y, Z)
   * ported from Proj 4.9.9 geocent.c
   */

  private HALF_PI = 0.5 * Math.PI;
  private COS_67P5 = 0.38268343236508977;
  private AD_C = 1.0026;

  private transform(srctype: number, dsttype: number, point: GeoPoint) {
    if (srctype === dsttype) return;

    if (
      (srctype !== 0 && srctype !== this.GRS80 && srctype !== this.UTMK) ||
      (dsttype !== 0 && dsttype !== this.GRS80 && dsttype !== this.UTMK)
    ) {
      this.geodetic_to_geocentric(srctype, point);

      if (srctype !== 0 && srctype !== this.GRS80 && srctype !== this.UTMK) {
        this.geocentric_to_wgs84(point);
      }

      if (dsttype !== 0 && dsttype !== this.GRS80 && dsttype !== this.UTMK) {
        this.geocentric_from_wgs84(point);
      }

      // Convert back to geodetic coordinates
      this.geocentric_to_geodetic(dsttype, point);
    }
  }

  private geodetic_to_geocentric(type: number, p: GeoPoint) {
    /*
     * The function Convert_Geodetic_To_Geocentric converts geodetic coordinates
     * (latitude, longitude, and height) to geocentric coordinates (X, Y, Z),
     * according to the current ellipsoid parameters.
     *
     *    Latitude  : Geodetic latitude in radians                     (input)
     *    Longitude : Geodetic longitude in radians                    (input)
     *    Height    : Geodetic height, in meters                       (input)
     *    X         : Calculated Geocentric X coordinate, in meters    (output)
     *    Y         : Calculated Geocentric Y coordinate, in meters    (output)
     *    Z         : Calculated Geocentric Z coordinate, in meters    (output)
     *
     */

    let Longitude = p.x;
    let Latitude = p.y;
    const Height = p.z;
    let X: number;
    let Y: number;
    let Z: number;

    let Rn: number;
    let SinLat: number;
    let Sin2Lat: number;
    let CosLat: number;

    /*
     ** Don't blow up if Latitude is just a little out of the value
     ** range as it may just be a rounding issue.  Also removed longitude
     ** test, it should be wrapped by Math.cos() and Math.sin().  NFW for PROJ.4, Sep/2001.
     */
    if (Latitude < -this.HALF_PI && Latitude > -1.001 * this.HALF_PI) {
      Latitude = -this.HALF_PI;
    } else if (Latitude > this.HALF_PI && Latitude < 1.001 * this.HALF_PI) {
      Latitude = this.HALF_PI;
    } else if (Latitude < -this.HALF_PI || Latitude > this.HALF_PI) {
      return true;
    }

    /* no erros */
    if (Longitude > Math.PI) {
      Longitude -= 2 * Math.PI;
    }

    SinLat = Math.sin(Latitude);
    CosLat = Math.cos(Latitude);
    Sin2Lat = SinLat * SinLat;
    Rn = this.mArMajor[type] / Math.sqrt(1.0 - this.m_Es[type] * Sin2Lat);
    X = (Rn + Height) * CosLat * Math.cos(Longitude);
    Y = (Rn + Height) * CosLat * Math.sin(Longitude);
    Z = (Rn * (1 - this.m_Es[type]) + Height) * SinLat;

    p.x = X;
    p.y = Y;
    p.z = Z;

    return false;
  }

  private geocentric_to_geodetic(type: number, p: GeoPoint) {
    const X = p.x;
    const Y = p.y;
    const Z = p.z;
    let Longitude: number;
    let Latitude = 0;
    let Height: number;

    let W: number; /* square of distance from Z axis */
    let W2: number; /* distance from Z axis */
    let T0: number; /* initial estimate of vertical component */
    let T1: number; /* corrected estimate of vertical component */
    let S0: number; /* initial estimate of horizontal component */
    let S1: number; /* corrected estimate of horizontal component */
    let SinB0: number; /* Math.sin(B0), B0 is estimate of Bowring aux doubleiable */
    let Sin3B0: number; /* cube of Math.sin(B0) */
    let CosB0: number; /* Math.cos(B0) */
    let SinP1: number; /* Math.sin(phi1), phi1 is estimated latitude */
    let CosP1: number; /* Math.cos(phi1) */
    let Rn: number; /* Earth radius at location */
    let Sum: number; /* numerator of Math.cos(phi1) */
    let AtPole: boolean; /* indicates location is in polar region */

    AtPole = false;

    if (X !== 0.0) {
      Longitude = Math.atan2(Y, X);
    } else {
      if (Y > 0) {
        Longitude = this.HALF_PI;
      } else if (Y < 0) {
        Longitude = -this.HALF_PI;
      } else {
        AtPole = true;
        Longitude = 0.0;

        if (Z > 0.0) {
          /* north pole */
          Latitude = this.HALF_PI;
        } else if (Z < 0.0) {
          /* south pole */
          Latitude = -this.HALF_PI;
        } else {
          /* center of earth */
          Latitude = this.HALF_PI;
          Height = -this.mArMinor[type];
          return;
        }
      }
    }

    W2 = X * X + Y * Y;
    W = Math.sqrt(W2);
    T0 = Z * this.AD_C;
    S0 = Math.sqrt(T0 * T0 + W2);
    SinB0 = T0 / S0;
    CosB0 = W / S0;
    Sin3B0 = SinB0 * SinB0 * SinB0;
    T1 = Z + this.mArMajor[type] * this.m_Esp[type] * Sin3B0;
    Sum = W - this.mArMajor[type] * this.m_Es[type] * CosB0 * CosB0 * CosB0;
    S1 = Math.sqrt(T1 * T1 + Sum * Sum);
    SinP1 = T1 / S1;
    CosP1 = Sum / S1;
    Rn = this.mArMajor[type] / Math.sqrt(1.0 - this.m_Es[type] * SinP1 * SinP1);

    if (CosP1 >= this.COS_67P5) {
      Height = W / CosP1 - Rn;
    } else if (CosP1 <= -this.COS_67P5) {
      Height = W / -CosP1 - Rn;
    } else {
      Height = Z / SinP1 + Rn * (this.m_Es[type] - 1.0);
    }

    if (!AtPole) {
      Latitude = Math.atan(SinP1 / CosP1);
    }

    p.x = Longitude;
    p.y = Latitude;
    p.z = Height;

    return;
  }

  private geocentric_to_wgs84(p: GeoPoint) {
    p.x += this.datumParams[0];
    p.y += this.datumParams[1];
    p.z += this.datumParams[2];
  }

  private geocentric_from_wgs84(p: GeoPoint) {
    p.x -= this.datumParams[0];
    p.y -= this.datumParams[1];
    p.z -= this.datumParams[2];
  }
}
