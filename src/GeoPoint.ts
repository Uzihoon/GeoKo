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

  set x(x: number) {
    this._x = x;
  }

  set y(y: number) {
    this._y = y;
  }

  set z(z: number) {
    this._z = z;
  }
}
