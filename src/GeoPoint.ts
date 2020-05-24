export default class GeoPoint {
  private _x: number;
  private _y: number;
  private _z: number;

  constructor(x?: number, y?: number, z?: number) {
    this._x = x || 0;
    this._y = y || 0;
    this._z = z || 0;
  }
  get x() {
    return this._x;
  }

  get y() {
    return this._y;
  }

  get z() {
    return this._z;
  }

  // tslint:disable-next-line: adjacent-overload-signatures
  set x(x: number) {
    this._x = x;
  }

  // tslint:disable-next-line: adjacent-overload-signatures
  set y(y: number) {
    this._y = y;
  }

  // tslint:disable-next-line: adjacent-overload-signatures
  set z(z: number) {
    this._z = z;
  }
}
