/* Logic/integration test of the offline controls. No WebGL/browser claim.
   Run: node scripts/check_viewer.cjs results/research_xcat/mapa_3d_interactivo.html */
const fs = require('fs'), vm = require('vm'), assert = require('assert/strict');
const html = fs.readFileSync(process.argv[2], 'utf8');
const payload = html.match(/<script type="application\/json" id="viewer-data">([\s\S]*?)<\/script>/)[1];
const data = JSON.parse(payload);
const code = html.slice(html.lastIndexOf('<script>')+8, html.lastIndexOf('</script>'));
new vm.Script(code); // syntax gate
class Element {
  constructor(id='') { this.id=id; this.value=''; this.checked=false; this.children=[]; this.handlers={}; this.style={}; this.dataset={}; this.textContent=''; }
  appendChild(child) { this.children.push(child); if(this.value==='' && child.value!==undefined)this.value=child.value; }
  replaceChildren() { this.children=[]; this.value=''; }
  addEventListener(event,fn) { this.handlers[event]=fn; }
  on(event,fn) { this.handlers[event]=fn; }
  fire(event) { if(this.handlers[event]) this.handlers[event]({target:this}); }
}
const elements = Object.fromEntries([...html.matchAll(/id="([^" ]+)"/g)].map(m=>[m[1],new Element(m[1])]));
elements['viewer-data'].textContent=payload;
elements['dose-map']._fullLayout={}; elements['point-index'].value='0'; elements['area-threshold'].value='20';
['x','y','z'].forEach(a=>elements['side-'+a].value='plus');
const views=['posterior','anterior','lateral','superior'].map(v=>{const b=new Element();b.dataset.view=v;return b;});
const traces=Array.from({length:13},()=>({})); const layouts=[]; let dah=null;
const Plotly={
  restyle(g,update,indices) { assert.equal(g,elements['dose-map']); indices.forEach(i=>{assert(i>=0&&i<13); for(const [key,value]of Object.entries(update))traces[i][key]=Array.isArray(value)?value[0]:value;}); },
  relayout(g,update) {layouts.push(update);},
  newPlot(id,series) {assert.equal(id,'dah');dah=series;}
};
const context={document:{getElementById(id){assert(elements[id],id);return elements[id];},createElement(){return new Element();},querySelectorAll(){return views;}},
  Plotly,window:{},setTimeout(fn){fn();return 1;},clearTimeout(){},console};
vm.runInNewContext(code,context);
assert(context.window.doseViewerReady); assert(dah);
const set=(id,value,event='change')=>{if(typeof value==='boolean')elements[id].checked=value;else elements[id].value=String(value);elements[id].fire(event);};
if(data.faces) assert.equal(traces[0].i.length,data.faces.length);
set('mode',1); set('show-beam',true); assert(traces[5].visible);
assert.deepEqual(Array.from(traces[6].x),data.geometry[0][0].source_patient_mm.slice(0,1).concat(
  data.geometry[0][0].isocenter_patient_mm[0],data.geometry[0][0].detector_patient_mm[0]));
set('show-gantry',true);assert.equal(traces[12].x.length,61);
set('inspect-points',true);assert(traces[8].customdata.length>0);
const point=traces[8].customdata[0];
elements['dose-map'].handlers.plotly_click({points:[{curveNumber:8,customdata:point}]});
assert(elements['point-info'].textContent.includes('Distancia a fuente'));
assert.equal(Number(elements['point-index'].value),point);
set('smooth',true);
if(data.faces){assert.equal(traces[0].facecolor.length,0);assert.equal(traces[0].vertexcolor.length,data.vertices.length);}
set('cut-z',0,'input');assert(traces[4].visible);
if(data.faces){assert(traces[0].z.every(z=>z>=-1e-6));assert.equal(traces[0].vertexcolor.length,traces[0].x.length);}
set('side-z','minus');if(data.faces)assert(traces[0].z.every(z=>z<=1e-6));
set('only-dose',true);set('mode',2);assert(elements.state.textContent.includes('Máximo'));
set('show-hotspots',true);assert(traces[11].visible);
set('area-threshold',20,'input');const expected=data.doses[0].reduce((sum,v,i)=>sum+(v>=.02?data.area_mm2[i]/100:0),0).toFixed(3);
assert(elements['area-result'].textContent.startsWith(expected));
elements['focus-dose'].fire('click');assert(layouts.at(-1)['scene.xaxis.range']);
elements['reset'].fire('click');assert(!elements['enable-z'].checked);assert(!elements.smooth.checked);
if(data.faces)assert.equal(traces[0].i.length,data.faces.length);
console.log(JSON.stringify({passed:true,scope:'DOM/Plotly stub: controls, clipping, event data, inspection, histogram; not WebGL rendering',
  face_count:data.faces?data.faces.length:null,point_count:data.centers_mm.length,inspected_point:point},null,2));
