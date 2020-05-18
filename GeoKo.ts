import GeoPoint from "./GeoPoint";

export default class GeoKo {
  public GEO = 0;
  public KATEC = 1;
  public TM = 2;
  public GRS80 = 3;
  public UTMK = 4;

  private m_Ind: number[] = [];
  private m_Es: number[] = [];
  private m_Esp: number[] = [];
  private src_m: number[] = [];
  private dst_m: number[] = [];

  private EPSLN = 0.0000000001;
  private m_arMajor: number[] = [];
  private m_arMinor: number[] = [];

  private m_arScaleFactor: number[] = [];
  private m_arLonCenter: number[] = [];
  private m_arLatCenter: number[] = [];
  private m_arFalseNorthing: number[] = [];
  private m_arFalseEasting: number[] = [];

  private datum_params: number[] = [];

  constructor() {
    // GEO Initial
    this.m_arScaleFactor[this.GEO] = 1;
    this.m_arLonCenter[this.GEO] = 0.0;
    this.m_arLatCenter[this.GEO] = 0.0;
    this.m_arFalseNorthing[this.GEO] = 0.0;
    this.m_arFalseEasting[this.GEO] = 0.0;
    this.m_arMajor[this.GEO] = 6378137.0;
    this.m_arMinor[this.GEO] = 6356752.3142;

    // KATEC Initial
    this.m_arScaleFactor[this.KATEC] = 0.9999;
    this.m_arLonCenter[this.KATEC] = 2.23402144255274; // 128
    this.m_arLatCenter[this.KATEC] = 0.663225115757845;
    this.m_arFalseNorthing[this.KATEC] = 600000.0;
    this.m_arFalseEasting[this.KATEC] = 400000.0;
    this.m_arMajor[this.KATEC] = 6377397.155;
    this.m_arMinor[this.KATEC] = 6356078.9633422494;

    // TM Initial
    this.m_arScaleFactor[this.TM] = 1.0;
    this.m_arLonCenter[this.TM] = 2.21661859489671;
    this.m_arLatCenter[this.TM] = 0.663225115757845;
    this.m_arFalseNorthing[this.TM] = 500000.0;
    this.m_arFalseEasting[this.TM] = 200000.0;
    this.m_arMajor[this.TM] = 6377397.155;
    this.m_arMinor[this.TM] = 6356078.9633422494;

    // GRS80 Initial
    this.m_arScaleFactor[this.GRS80] = 1.0;
    this.m_arLonCenter[this.GRS80] = 2.2165681500328;
    this.m_arLatCenter[this.GRS80] = 0.663225115757845;
    this.m_arFalseNorthing[this.GRS80] = 500000.0;
    this.m_arFalseEasting[this.GRS80] = 200000.0;
    this.m_arMajor[this.GRS80] = 6378137;
    this.m_arMinor[this.GRS80] = 6356752.3142;

    // UTMK Initial
    this.m_arScaleFactor[this.UTMK] = 0.9996;
    this.m_arLonCenter[this.UTMK] = 2.22529479629277;
    this.m_arLatCenter[this.UTMK] = 0.663225115757845;
    this.m_arFalseNorthing[this.UTMK] = 2000000.0;
    this.m_arFalseEasting[this.UTMK] = 1000000.0;
    this.m_arMajor[this.UTMK] = 6378137;
    this.m_arMinor[this.UTMK] = 6356752.3141403558;

    this.datum_params[0] = -146.43;
    this.datum_params[1] = 507.89;
    this.datum_params[2] = 681.46;

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
    let tmp = this.m_arMinor[type] / this.m_arMajor[type];

    this.m_Es[type] = 1.0 - tmp * tmp;
    this.m_Esp[type] = this.m_Es[type] / (1.0 - this.m_Es[type]);

    if (this.m_Es[type] < 0.00001) {
      this.m_Ind[type] = 1.0;
    } else {
      this.m_Ind[type] = 0.0;
    }
  }

  private initialSrcDst(type: number) {
    const es = this.m_Es[type];
    const data =
      this.m_arMajor[type] *
      this.mlfn(
        this.e0fn(es),
        this.e1fn(es),
        this.e2fn(es),
        this.e3fn(es),
        this.m_arLatCenter[type]
      );
    this.src_m[type] = data;
    this.dst_m[type] = data;
  }

  private D2R(degree: number) {
    return (degree * Math.PI) / 180.0;
  }

  private R2D(radian: number) {
    return (radian * 180.0) / Math.PI;
  }

  private mlfn(e0: number, e1: number, e2: number, e3: number, phi: number) {
    return (
      e0 * phi -
      e1 * Math.sin(2.0 * phi) +
      e2 * Math.sin(4.0 * phi) -
      e3 * Math.sin(6.0 * phi)
    );
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

  public convert(srctype: number, dsttype: number, in_pt: GeoPoint) {
    const tmpPt = new GeoPoint();
    const out_pt = new GeoPoint();

    if (srctype === this.GEO) {
      tmpPt.x = this.D2R(in_pt.x);
      tmpPt.y = this.D2R(in_pt.y);
    } else {
      this.tm2geo(srctype, in_pt, tmpPt);
    }

    if (dsttype === this.GEO) {
      out_pt.x = this.R2D(tmpPt.x);
      out_pt.y = this.R2D(tmpPt.y);
    } else {
      this.geo2tm(dsttype, tmpPt, out_pt);
    }

    return out_pt;
  }

  public geo2tm(dsttype: number, in_pt: GeoPoint, out_pt: GeoPoint) {
    let x: number;
    let y: number;

    this.transform(this.GEO, dsttype, in_pt);

    let delta_lon = in_pt.x - this.m_arLonCenter[dsttype];
    let sin_phi = Math.sin(in_pt.y);
    let cos_phi = Math.cos(in_pt.y);

    if (this.m_Ind[dsttype] !== 0) {
      let b = cos_phi * Math.sin(delta_lon);

      if (Math.abs(Math.abs(b) - 1.0) < this.EPSLN) {
        // 무한대 에러
        console.error("무한대 에러");
      }
    } else {
      let b = 0;
      x =
        0.5 *
        this.m_arMajor[dsttype] *
        this.m_arScaleFactor[dsttype] *
        Math.log((1.0 + b) / (1.0 - b));

      let con = Math.acos(
        (cos_phi * Math.cos(delta_lon)) / Math.sqrt(1.0 - b * b)
      );

      if (in_pt.y < 0) {
        con = con * -1;
        y =
          this.m_arMajor[dsttype] *
          this.m_arScaleFactor[dsttype] *
          (con - this.m_arLatCenter[dsttype]);
      }
    }

    const al = cos_phi * delta_lon;
    const als = al * al;
    const c = this.m_Es[dsttype] * cos_phi * cos_phi;
    const tq = Math.tan(in_pt.y);
    const t = tq * tq;
    const con = 1.0 - this.m_Es[dsttype] * sin_phi * sin_phi;
    const n = this.m_arMajor[dsttype] / Math.sqrt(con);
    const ml =
      this.m_arMajor[dsttype] *
      this.mlfn(
        this.e0fn(this.m_Es[dsttype]),
        this.e1fn(this.m_Es[dsttype]),
        this.e2fn(this.m_Es[dsttype]),
        this.e3fn(this.m_Es[dsttype]),
        in_pt.y
      );

    out_pt.x =
      this.m_arScaleFactor[dsttype] *
        n *
        al *
        (1.0 +
          (als / 6.0) *
            (1.0 -
              t +
              c +
              (als / 20.0) *
                (5.0 -
                  18.0 * t +
                  t * t +
                  72.0 * c -
                  58.0 * this.m_Esp[dsttype]))) +
      this.m_arFalseEasting[dsttype];

    out_pt.y =
      this.m_arScaleFactor[dsttype] *
        (ml -
          this.dst_m[dsttype] +
          n *
            tq *
            (als *
              (0.5 +
                (als / 24.0) *
                  (5.0 -
                    t +
                    9.0 * c +
                    4.0 * c * c +
                    (als / 30.0) *
                      (61.0 -
                        58.0 * t +
                        t * t +
                        600.0 * c -
                        330.0 * this.m_Esp[dsttype]))))) +
      this.m_arFalseNorthing[dsttype];
  }

  public tm2geo(srctype: number, in_pt: GeoPoint, out_pt: GeoPoint) {
    const tmpPt = new GeoPoint(in_pt.x, in_pt.y);
    const max_iter = 6;

    if (this.m_Ind[srctype] !== 0) {
      const f = Math.exp(
        in_pt.x / (this.m_arMajor[srctype] * this.m_arScaleFactor[srctype])
      );
      const g = 0.5 * (f - 1.0 / f);
      const temp =
        this.m_arLatCenter[srctype] +
        tmpPt.y / (this.m_arMajor[srctype] * this.m_arScaleFactor[srctype]);
      const h = Math.cos(temp);
      const con = Math.sqrt((1.0 - h * h) / (1.0 + g * g));
      out_pt.y = this.asinz(con);

      if (temp < 0) out_pt.y *= -1;

      if (g === 0 && h === 0) {
        out_pt.x = this.m_arLonCenter[srctype];
      } else {
        out_pt.x = Math.atan(g / h) + this.m_arLonCenter[srctype];
      }
    }

    tmpPt.x -= this.m_arFalseEasting[srctype];
    tmpPt.y -= this.m_arFalseNorthing[srctype];

    const con =
      (this.src_m[srctype] + tmpPt.y / this.m_arScaleFactor[srctype]) /
      this.m_arMajor[srctype];

    let phi = con;
    let i = 0;

    while (true) {
      const delta_Phi =
        (con +
          this.e1fn(this.m_Es[srctype]) * Math.sin(2.0 * phi) -
          this.e2fn(this.m_Es[srctype]) * Math.sin(4.0 * phi) +
          this.e3fn(this.m_Es[srctype]) * Math.sin(6.0 * phi)) /
          this.e0fn(this.m_Es[srctype]) -
        phi;

      phi = phi + delta_Phi;

      if (Math.abs(delta_Phi) <= this.EPSLN) break;

      if (i >= max_iter) {
        break;
      }

      i++;
    }

    if (Math.abs(phi) < Math.PI / 2) {
      const sin_phi = Math.sin(phi);
      const cos_phi = Math.cos(phi);
      const tan_phi = Math.tan(phi);
      const c = this.m_Esp[srctype] * cos_phi * cos_phi;
      const cs = c * c;
      const t = tan_phi * tan_phi;
      const ts = t * t;
      const cont = 1.0 - this.m_Es[srctype] * sin_phi * sin_phi;
      const n = this.m_arMajor[srctype] / Math.sqrt(cont);
      const r = (n * (1.0 - this.m_Es[srctype])) / cont;
      const d = tmpPt.x / (n * this.m_arScaleFactor[srctype]);
      const ds = d * d;

      out_pt.y =
        phi -
        ((n * tan_phi * ds) / r) *
          (0.5 -
            (ds / 24.0) *
              (5.0 +
                3.0 * t +
                10.0 * c -
                4.0 * cs -
                9.0 * this.m_Esp[srctype] -
                (ds / 30.0) *
                  (61.0 +
                    90.0 * t +
                    298.0 * c +
                    45.0 * ts -
                    252.0 * this.m_Esp[srctype] -
                    3.0 * cs)));
      out_pt.x =
        this.m_arLonCenter[srctype] +
        (d *
          (1.0 -
            (ds / 6.0) *
              (1.0 +
                2.0 * t +
                c -
                (ds / 20.0) *
                  (5.0 -
                    2.0 * c +
                    28.0 * t -
                    3.0 * cs +
                    8.0 * this.m_Esp[srctype] +
                    24.0 * ts)))) /
          cos_phi;
    } else {
      out_pt.y = Math.PI * 0.5 * Math.sin(tmpPt.y);
      out_pt.x = this.m_arLonCenter[srctype];
    }

    this.transform(srctype, this.GEO, out_pt);
  }

  public getDistancebyGeo(pt1: GeoPoint, pt2: GeoPoint) {
    const lat1 = this.D2R(pt1.y);
    const lon1 = this.D2R(pt1.x);
    const lat2 = this.D2R(pt2.y);
    const lon2 = this.D2R(pt2.x);

    const longitude = lon2 - lon1;
    const latitude = lat2 - lat1;

    const a =
      Math.pow(Math.sin(latitude / 2.0), 2) +
      Math.cos(lat1) * Math.cos(lat2) * Math.pow(Math.sin(longitude / 2.0), 2);

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
    let Height = p.z;
    let X: number;
    let Y: number;
    let Z: number;

    let Rn: number;
    let Sin_Lat: number;
    let Sin2_Lat: number;
    let Cos_Lat: number;

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

    Sin_Lat = Math.sin(Latitude);
    Cos_Lat = Math.cos(Latitude);
    Sin2_Lat = Sin_Lat * Sin_Lat;
    Rn = this.m_arMajor[type] / Math.sqrt(1.0 - this.m_Es[type] * Sin2_Lat);
    X = (Rn + Height) * Cos_Lat * Math.cos(Longitude);
    Y = (Rn + Height) * Cos_Lat * Math.sin(Longitude);
    Z = (Rn * (1 - this.m_Es[type]) + Height) * Sin_Lat;

    p.x = X;
    p.y = Y;
    p.z = Z;

    return false;
  }

  private geocentric_to_geodetic(type: number, p: GeoPoint) {
    let X = p.x;
    let Y = p.y;
    let Z = p.z;
    let Longitude: number;
    let Latitude = 0;
    let Height: number;

    let W: number; /* square of distance from Z axis */
    let W2: number; /* distance from Z axis */
    let T0: number; /* initial estimate of vertical component */
    let T1: number; /* corrected estimate of vertical component */
    let S0: number; /* initial estimate of horizontal component */
    let S1: number; /* corrected estimate of horizontal component */
    let Sin_B0: number; /* Math.sin(B0), B0 is estimate of Bowring aux doubleiable */
    let Sin3_B0: number; /* cube of Math.sin(B0) */
    let Cos_B0: number; /* Math.cos(B0) */
    let Sin_p1: number; /* Math.sin(phi1), phi1 is estimated latitude */
    let Cos_p1: number; /* Math.cos(phi1) */
    let Rn: number; /* Earth radius at location */
    let Sum: number; /* numerator of Math.cos(phi1) */
    let At_Pole: boolean; /* indicates location is in polar region */

    At_Pole = false;

    if (X != 0.0) {
      Longitude = Math.atan2(Y, X);
    } else {
      if (Y > 0) {
        Longitude = this.HALF_PI;
      } else if (Y < 0) {
        Longitude = -this.HALF_PI;
      } else {
        At_Pole = true;
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
          Height = -this.m_arMinor[type];
          return;
        }
      }
    }

    W2 = X * X + Y * Y;
    W = Math.sqrt(W2);
    T0 = Z * this.AD_C;
    S0 = Math.sqrt(T0 * T0 + W2);
    Sin_B0 = T0 / S0;
    Cos_B0 = W / S0;
    Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;
    T1 = Z + this.m_arMajor[type] * this.m_Esp[type] * Sin3_B0;
    Sum = W - this.m_arMajor[type] * this.m_Es[type] * Cos_B0 * Cos_B0 * Cos_B0;
    S1 = Math.sqrt(T1 * T1 + Sum * Sum);
    Sin_p1 = T1 / S1;
    Cos_p1 = Sum / S1;
    Rn =
      this.m_arMajor[type] / Math.sqrt(1.0 - this.m_Es[type] * Sin_p1 * Sin_p1);

    if (Cos_p1 >= this.COS_67P5) {
      Height = W / Cos_p1 - Rn;
    } else if (Cos_p1 <= -this.COS_67P5) {
      Height = W / -Cos_p1 - Rn;
    } else {
      Height = Z / Sin_p1 + Rn * (this.m_Es[type] - 1.0);
    }

    if (!At_Pole) {
      Latitude = Math.atan(Sin_p1 / Cos_p1);
    }

    p.x = Longitude;
    p.y = Latitude;
    p.z = Height;

    return;
  }

  private geocentric_to_wgs84(p: GeoPoint) {
    p.x += this.datum_params[0];
    p.y += this.datum_params[1];
    p.z += this.datum_params[2];
  }

  private geocentric_from_wgs84(p: GeoPoint) {
    p.x -= this.datum_params[0];
    p.y -= this.datum_params[1];
    p.z -= this.datum_params[2];
  }
}
