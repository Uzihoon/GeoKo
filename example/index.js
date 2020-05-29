/**
 * handle 3D rendering
 */

(function renderingEarth() {
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
  const renderer = new THREE.WebGLRenderer();

  renderer.setSize(window.innerWidth, window.innerHeight);
  document.body.appendChild(renderer.domElement);

  /* Create Lights */
  const spotLight = new THREE.SpotLight(0xffffff);
  spotLight.position.set(100, 100, 100);
  spotLight.castShadow = true;
  spotLight.shadowMapWidth = 1024;
  spotLight.shadowMapHeight = 1024;
  spotLight.shadowCameraNear = 500;
  spotLight.shadowCameraFar = 4000;
  spotLight.shadowCameraFov = 30;
  scene.add(spotLight);

  /* Create Material */
  function Mat() {
    const material = new THREE.MeshPhongMaterial({
      color: new THREE.Color('rgb(0, 76, 255)'),
      emissive: new THREE.Color('rgb(0, 117, 255)'),
      specular: new THREE.Color('rgb(0, 141, 255)'),
      shininess: 1,
      shading: THREE.FlatShading,
      wireframe: 1,
      transparent: 1,
      opacity: 0.15,
    });

    return material;
  }

  /* Create Geometry */
  const geometry = new THREE.SphereGeometry(50, 20, 20, 0, Math.PI * 2, 0, Math.PI);

  /* Create Earth Sphere */
  const earth = new THREE.Mesh(geometry, Mat());

  scene.add(earth);

  camera.position.z = 90;

  function render() {
    requestAnimationFrame(render);
    earth.rotation.x += 0.001;
    earth.rotation.y += 0.001;
    renderer.render(scene, camera);
  }

  render();
})();

/**
 * Geolocation
 */

(function handleGeoKo() {
  const getCode = (el) => el.dataset.convert;
  const setGeo = (id, code) => (document.getElementById(id).innerHTML = code);
  const getGeo = (id) => +document.getElementById(id).value || 0;
  const remove = (t) => t.removeChild(t.children[0]);

  function handleRemove(e) {
    e.preventDefault();
    e.stopPropagation();

    const fromEl = document.getElementById('from');
    const toEl = document.getElementById('to');

    if (fromEl.children.length) remove(fromEl);
    if (toEl.children.length) remove(toEl);

    window.removeEventListener('click', handleRemove);
  }

  /**
   * Handle Select Box
   */
  (function handleSelectBox() {
    const convertList = [
      {
        id: 'geo',
        label: 'GEO',
      },
      {
        id: 'tm',
        label: 'TM',
      },
      {
        id: 'katec',
        label: 'KATEC',
      },
      {
        id: 'grs80',
        label: 'GRS80',
      },
      {
        id: 'utmk',
        label: 'UTMK',
      },
    ];

    const fromEl = document.getElementById('from');
    const toEl = document.getElementById('to');

    const onClickConvert = (e, t, id, label) => {
      e.preventDefault();
      e.stopPropagation();
      t.dataset.convert = id;
      t.innerText = label;
    };

    function onClick(e) {
      e.preventDefault();
      e.stopPropagation();

      const code = getCode(this);
      const that = this;

      if (this.children.length) {
        remove(this);
        return;
      }
      const select = document.createElement('div');
      select.className = 'select-list';

      convertList.map((convert) => {
        if (convert.id !== code) {
          const convertEl = document.createElement('div');
          convertEl.className = 'convert-list';
          convertEl.id = convert.id;
          convertEl.innerHTML = convert.label;
          convertEl.addEventListener('click', (e) => onClickConvert(e, that, convert.id, convert.label));
          select.appendChild(convertEl);
        }
      });

      this.appendChild(select);

      window.addEventListener('click', handleRemove);
    }

    fromEl.addEventListener('click', onClick);
    toEl.addEventListener('click', onClick);
  })();

  /**
   * Handle Convert
   */
  (function handleConvert() {
    const convertBtn = document.getElementById('convert');

    const convertGeolocation = () => {
      const from = document.getElementById('from');
      const to = document.getElementById('to');

      if (getCode(from) === getCode(to)) return;

      const geoko = new GeoKo();

      const converted = geoko.setX(getGeo('x_axis')).setY(getGeo('y_axis'))[`${getCode(from)}_to_${getCode(to)}`]();

      setGeo('converted_from', converted.x);
      setGeo('converted_to', converted.y);
    };

    convertBtn.addEventListener('click', convertGeolocation);
  })();
})();
