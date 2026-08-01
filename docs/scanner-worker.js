importScripts('cmbs-decoder.js');

var CMBS = self.CMBS;

self.onmessage = function (e) {
  var msg = e.data;
  if (msg.type !== 'decode') return;
  var t0 = self.performance.now();
  var rgb = new Uint8Array(msg.width * msg.height * 3);
  var rgba = msg.buf;
  for (var i = 0, o = 0; i < msg.width * msg.height; i++, o += 3) {
    var j = i * 4;
    rgb[o] = rgba[j];
    rgb[o + 1] = rgba[j + 1];
    rgb[o + 2] = rgba[j + 2];
  }
  var result = null;
  try {
    result = CMBS.decodeFrame(rgb, msg.width, msg.height);
  } catch (err) {
    result = { error: String(err && err.message || err) };
  }
  self.postMessage({
    type: 'result',
    id: msg.id,
    result: result,
    ms: self.performance.now() - t0,
  }, []);
};
