import GeoPoint from './GeoPoint';
import ConvertLocation from './ConvertLocation';

export default class GeoKo {
  private point: GeoPoint;
  private converter: ConvertLocation;

  constructor() {
    this.point = new GeoPoint(0, 0);
    this.converter = new ConvertLocation();
  }

  public setX(x: number) {
    this.point.x = x;
    return this;
  }

  public setY(y: number) {
    this.point.y = y;
    return this;
  }

  public utmk_to_geo() {
    return this.convert(4, 0);
  }

  public utmk_to_tm() {
    return this.convert(4, 2);
  }

  public utmk_to_katec() {
    return this.convert(4, 1);
  }

  public utmk_to_grs80() {
    return this.convert(4, 3);
  }

  public geo_to_tm() {
    return this.convert(0, 2);
  }

  public geo_to_katec() {
    return this.convert(0, 1);
  }

  public geo_to_grs80() {
    return this.convert(0, 3);
  }

  public geo_to_utmk() {
    return this.convert(0, 4);
  }

  public grs80_to_geo() {
    return this.convert(3, 0);
  }

  public grs80_to_tm() {
    return this.convert(3, 2);
  }

  public grs80_to_katec() {
    return this.convert(3, 1);
  }

  public grs80_to_utmk() {
    return this.convert(3, 4);
  }

  public tm_to_geo() {
    return this.convert(2, 0);
  }

  public tm_to_katec() {
    return this.convert(2, 1);
  }

  public tm_to_grs80() {
    return this.convert(2, 3);
  }

  public tm_to_utmk() {
    return this.convert(2, 4);
  }

  public katec_to_geo() {
    return this.convert(1, 0);
  }

  public katec_to_tm() {
    return this.convert(1, 2);
  }

  public katec_to_grs80() {
    return this.convert(1, 3);
  }

  public katec_to_utmk() {
    return this.convert(1, 4);
  }

  private convert(geoType: number, convertType: number) {
    return this.converter.convert(geoType, convertType, this.point);
  }
}
