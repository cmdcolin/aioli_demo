(this.webpackJsonpaioli_demo=this.webpackJsonpaioli_demo||[]).push([[0],{20:function(t,e,i){},37:function(t,e,i){"use strict";i.r(e);var s=i(0),a=i.n(s),n=i(15),c=i.n(n),r=(i(20),i(3)),o=i.n(r),l=i(7),h=i(2),f=i(6),u=i(10),b=i(11),d=1e4,m=2e5,p=function(){function t(e){Object(u.a)(this,t),this.rowNumber=e,this.padding=1,this.sizeLimit=1e7}return Object(b.a)(t,[{key:"log",value:function(t){console.log("r".concat(this.rowNumber,": ").concat(t))}},{key:"setAllFilled",value:function(t){this.allFilled=t}},{key:"getItemAt",value:function(t){if(this.allFilled)return this.allFilled;if(void 0!==this.min&&!(t<this.min)&&!(t>=this.max)){var e=t-this.offset;return this.bits[e]}}},{key:"isRangeClear",value:function(t,e){if(this.allFilled)return!1;if(void 0===this.min)return!0;if(e<=this.min||t>=this.max)return!0;for(var i=Math.min(this.max,e),s=Math.max(this.min,t);s<e&&s<i;s+=1)if(this.getItemAt(s))return!1;return!0}},{key:"initialize",value:function(t,e){var i=e-t;this.offset=t-i,this.min=t,this.max=e,this.bits=new Array(e-t+2*i)}},{key:"addRect",value:function(t,e){var i=t.l,s=t.r+this.padding;if(void 0===this.min)this.initialize(i,s);else{var a=this.bits.length;if(s-this.offset>=this.bits.length){var n=s-this.offset-this.bits.length+1+this.bits.length;this.bits.length+n>this.sizeLimit?(console.warn("Layout width limit exceeded, discarding old layout. Please be more careful about discarding unused blocks."),this.initialize(i,s)):n>0&&(this.bits=this.bits.concat(new Array(n)))}if(i<this.offset){var c=this.offset-i+a;this.bits.length+c>this.sizeLimit?(console.warn("Layout width limit exceeded, discarding old layout. Please be more careful about discarding unused blocks."),this.initialize(i,s)):(this.bits=new Array(c).concat(this.bits),this.offset-=c)}}var r=i-this.offset,o=s-this.offset;o-r>m&&console.warn("Layout X pitch set too low, feature spans ".concat(o-r," bits in a single row."),t,e);for(var l=r;l<o;l+=1)this.bits[l]=e;i<this.min&&(this.min=i),s>this.max&&(this.max=s)}},{key:"discardRange",value:function(t,e){if(!this.allFilled&&this.bits&&!(e<=this.min||t>=this.max)){if(t<=this.min&&e>=this.max)return this.min=void 0,this.max=void 0,this.bits=void 0,void(this.offset=void 0);if(e>this.min&&t<=this.min&&(this.min=e),t<this.max&&e>=this.max&&(this.max=t),this.offset<this.min-d&&this.bits.length>this.max+d-this.offset){var i=this.min-this.offset,s=this.bits.length-1-(this.max-this.offset);this.bits=this.bits.slice(i,this.bits.length-s),this.offset+=i}else if(this.offset<this.min-d){var a=this.min-Math.floor(5e3)-this.offset;this.bits.splice(0,a),this.offset+=a}else if(this.bits.length>this.max-this.offset+d){var n=this.max-this.offset+1+Math.floor(5e3);this.bits.length=n}for(var c=Math.max(this.min,t)-this.offset,r=Math.min(e,this.max)-this.offset,o=c;o>=0&&o<r;o+=1)this.bits[o]=void 0}}}]),t}(),v=function(){function t(){var e=arguments.length>0&&void 0!==arguments[0]?arguments[0]:{};Object(u.a)(this,t),this.pitchX=e.pitchX||5,this.pitchY=e.pitchY||3,this.displayMode=e.displayMode,"compact"===this.displayMode&&(this.pitchY=Math.round(this.pitchY/4)||1,this.pitchX=Math.round(this.pitchX/4)||1),this.bitmap=[],this.rectangles={},this.maxHeight=Math.ceil((e.maxHeight||1/0)/this.pitchY),this.pTotalHeight=0}return Object(b.a)(t,[{key:"addRect",value:function(t,e,i,s,a){if(t in this.rectangles){var n=this.rectangles[t];return null===n.top?null:(this._addRectToBitmap(n,a),n.top*this.pitchY)}var c=Math.floor(e/this.pitchX),r=Math.floor(i/this.pitchX),o=Math.ceil(s/this.pitchY),l={id:t,l:c,r:r,mX:Math.floor((c+r)/2),h:o};a&&(l.data=a);for(var h=this.maxHeight-o,f=0;f<=h&&this._collides(l,f);f+=1);return f>h?(l.top=f=null,this.rectangles[t]=l,this.pTotalHeight=Math.max(this.pTotalHeight||0,f+o),null):(l.top=f,this._addRectToBitmap(l,a),this.rectangles[t]=l,this.pTotalHeight=Math.max(this.pTotalHeight||0,f+o),f*this.pitchY)}},{key:"_collides",value:function(t,e){if("collapsed"===this.displayMode)return!1;for(var i=this.bitmap,s=e+t.h,a=e;a<s;a+=1){var n=i[a];if(n&&!n.isRangeClear(t.l,t.r))return!0}return!1}},{key:"_autovivifyRow",value:function(t,e){var i=t[e];return i||(i=new p(e),t[e]=i),i}},{key:"_addRectToBitmap",value:function(t,e){if(null!==t.top){e=e||!0;var i=this.bitmap,s=this._autovivifyRow,a=t.top+t.h;if(t.r-t.l>m)for(var n=t.top;n<a;n+=1)s(i,n).setAllFilled(e);else for(var c=t.top;c<a;c+=1)s(i,c).addRect(t,e)}}},{key:"discardRange",value:function(t,e){for(var i=Math.floor(t/this.pitchX),s=Math.floor(e/this.pitchX),a=this.bitmap,n=0;n<a.length;n+=1){var c=a[n];c&&c.discardRange(i,s)}}},{key:"hasSeen",value:function(t){return!!this.rectangles[t]}},{key:"getByCoord",value:function(t,e){var i=Math.floor(e/this.pitchY),s=this.bitmap[i];if(s){var a=Math.floor(t/this.pitchX);return s.getItemAt(a)}}},{key:"getByID",value:function(t){var e=this.rectangles[t];if(e)return e.data||!0}},{key:"cleanup",value:function(){}},{key:"getTotalHeight",value:function(){return this.pTotalHeight*this.pitchY}}]),t}(),j=i(5),x=i(1);function g(t){var e=t.split(":"),i=Object(h.a)(e,2),s=i[0],a=i[1].split("-"),n=Object(h.a)(a,2);return{refId:+s-1,start:+n[0],end:+n[1]}}var O=1800,y=100;var w=function(){var t=Object(s.useRef)(),e=Object(s.useRef)(),i=Object(s.useRef)(),a=Object(j.c)({loc:Object(j.d)(j.b,"1:20000-40000"),file:Object(j.d)(j.b,"https://s3.amazonaws.com/1000genomes/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam")}),n=Object(h.a)(a,2),c=n[0],r=n[1],u=Object(s.useState)(),b=Object(h.a)(u,2),d=b[0],m=b[1],p=Object(s.useState)(),w=Object(h.a)(p,2),k=w[0],M=w[1],R=Object(s.useState)(),A=Object(h.a)(R,2),S=A[0],F=A[1],H=Object(s.useState)(),L=Object(h.a)(H,2),T=L[0],C=L[1],I=Object(s.useState)(c.file),X=Object(h.a)(I,2),B=X[0],E=X[1],Y=Object(s.useState)(c.loc),_=Object(h.a)(Y,2),z=_[0],N=_[1],P=Object(s.useState)(),W=Object(h.a)(P,2),D=W[0],G=W[1],U=function(){var t=Object(s.useState)(0),e=Object(h.a)(t,2)[1];return Object(s.useCallback)((function(){e((function(t){return t+1}))}),[])}(),J=Object(s.useState)(),q=Object(h.a)(J,2),K=q[0],Q=q[1];return Object(s.useEffect)((function(){Object(l.a)(o.a.mark((function t(){var e;return o.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:return e=new f.Aioli("samtools/1.10"),t.next=3,e.init();case 3:F(e);case 4:case"end":return t.stop()}}),t)})))()}),[]),Object(s.useEffect)((function(){Object(l.a)(o.a.mark((function t(){var e,i,s,a;return o.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:if(!S){t.next=24;break}if(!K){t.next=11;break}return e=K[0].name.endsWith("bam")||K[0].name.endsWith("cram")?0:1,t.next=5,f.Aioli.mount(K[e]);case 5:return i=t.sent,t.next=8,f.Aioli.mount(K[Number(!e)]);case 8:C(i),t.next=24;break;case 11:return s=new URL(c.file,window.location),t.next=14,f.Aioli.mount("".concat(s));case 14:if(a=t.sent,!c.file.endsWith("bam")){t.next=20;break}return t.next=18,f.Aioli.mount("".concat(s,".bai"));case 18:t.next=23;break;case 20:if(!c.file.endsWith("cram")){t.next=23;break}return t.next=23,f.Aioli.mount("".concat(s,".crai"));case 23:C(a);case 24:case"end":return t.stop()}}),t)})))()}),[c.file,S,c.fasta,K]),Object(s.useEffect)((function(){Object(l.a)(o.a.mark((function t(){var e;return o.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:if(!T||!S){t.next=5;break}return t.next=3,S.exec("view ".concat(T.path," ").concat(c.loc));case 3:e=t.sent,m(e);case 5:case"end":return t.stop()}}),t)})))()}),[T,c.loc,S]),Object(s.useEffect)((function(){Object(l.a)(o.a.mark((function t(){var e;return o.a.wrap((function(t){for(;;)switch(t.prev=t.next){case 0:if(!T||!S){t.next=5;break}return t.next=3,S.exec("mpileup -r ".concat(c.loc," ").concat(T.path));case 3:e=t.sent,M(e);case 5:case"end":return t.stop()}}),t)})))()}),[T,c.loc,S]),Object(s.useEffect)((function(){var e=t.current.getContext("2d");e.clearRect(0,0,O,1e3);var i=g(c.loc),s=O/(i.end-i.start),a=new v;d&&!d.stdout&&G(d.stderr),null===d||void 0===d||d.stdout.split("\n").filter((function(t){return!!t})).forEach((function(t,n){for(var c=t.split("\t"),r=Object(h.a)(c,6),o=r[1],l=r[3],f=r[5],u=+l-1,b=+o,d=(f||"").split(/([MIDNSHPX=])/),m=0,p=0;p<d.length;p+=2){var v=+d[p],j=d[p+1];"I"!==j&&"S"!==j&&"H"!==j&&(m+=v)}e.fillStyle=16&b?"#99f":"#f99";var x=u+m,g=(u-i.start)*s,O=(x-u)*s,y=a.addRect(n,u,x,10);e.fillRect(g,y,O,10)}))}),[d,c.loc]),Object(s.useEffect)((function(){var t=i.current.getContext("2d");t.clearRect(0,0,O,y);var e=g(c.loc),s=O/(e.end-e.start),a=0;k&&(k.stdout.split("\n").forEach((function(t){var e=t.split("\t"),i=+Object(h.a)(e,4)[3];a=Math.max(a,i||0)})),k.stdout.split("\n").forEach((function(i){var n=i.split("\t"),c=Object(h.a)(n,4),r=+c[1]-1,o=+c[3],l=(r-e.start)*s,f=(r+1-r)*s;t.fillStyle="#ccc";var u=o/a*y;t.fillRect(l,y-u,f+.9,u)})),t.fillStyle="black",t.fillText("[0, ".concat(a,"]"),0,20))}),[k,c.loc]),Object(x.jsxs)("div",{children:[Object(x.jsx)("p",{children:"Enter BAM/CRAM file and location. This app uses @biowasm/aioli's packaging of samtools to run samtools view and samtools mpileup"}),Object(x.jsxs)("form",{onSubmit:function(t){r({file:B,loc:z}),e.current.files.length&&Q(e.current.files),M(),m(),G(),U(),t.preventDefault()},children:[Object(x.jsx)("label",{htmlFor:"url",children:"URL: "}),Object(x.jsx)("input",{id:"url",type:"text",value:B,style:{minWidth:"75%"},onChange:function(t){return E(t.target.value)}}),Object(x.jsx)("br",{}),Object(x.jsx)("label",{htmlFor:"file",children:"File (import both BAM and BAI): "}),Object(x.jsx)("input",{id:"file",ref:e,type:"file",multiple:"multiple"}),Object(x.jsx)("br",{}),Object(x.jsx)("label",{htmlFor:"loc",children:"Location: "}),Object(x.jsx)("input",{id:"loc",type:"text",value:z,onChange:function(t){return N(t.target.value)}}),Object(x.jsx)("button",{type:"submit",children:"Submit"})]}),d?null:Object(x.jsx)("div",{className:"dots",children:"Loading..."}),D?Object(x.jsx)("div",{style:{color:"red"},children:D}):null,Object(x.jsx)("canvas",{ref:i,width:O,height:y}),Object(x.jsx)("canvas",{ref:t,width:O,height:1e3})]})},k=function(t){t&&t instanceof Function&&i.e(3).then(i.bind(null,38)).then((function(e){var i=e.getCLS,s=e.getFID,a=e.getFCP,n=e.getLCP,c=e.getTTFB;i(t),s(t),a(t),n(t),c(t)}))};c.a.render(Object(x.jsx)(a.a.StrictMode,{children:Object(x.jsx)(j.a,{children:Object(x.jsx)(w,{})})}),document.getElementById("root")),k()}},[[37,1,2]]]);
//# sourceMappingURL=main.49ed3ba0.chunk.js.map