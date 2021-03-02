import { useEffect, useRef, useState } from "react";
import { Aioli } from "@biowasm/aioli";

function parseLocString(locString) {
  const [refId, rest] = locString.split(":");
  const [start, end] = rest.split("-");
  return { refId: +refId - 1, start: +start, end: +end };
}

export function parseCigar(cigar: string) {
  return (cigar || "").split(/([MIDNSHPX=])/);
}

const featHeight = 10;
const spacing = 1;
const width = 1800;
const height = 1000;
const initialLoc = "ctgA:1-1000";
const initialFile = "volvox-sorted.bam";
function App() {
  const [data, setData] = useState();
  const ref = useRef();

  const [file, setFile] = useState(initialFile);
  const [loc, setLoc] = useState(initialLoc);
  const [fileVal, setFileVal] = useState(initialFile);
  const [locVal, setLocVal] = useState(initialLoc);
  const [samtools, setSamtools] = useState();

  useEffect(() => {
    (async () => {
      const tool = new Aioli("samtools/1.10");
      await tool.init();
      setSamtools(tool);
    })();
  }, []);

  useEffect(() => {
    (async () => {
      if (samtools) {
        const url = new URL(fileVal, window.location);
        const bam = await Aioli.mount(`${url}`);
        const bai = await Aioli.mount(`${url}.bai`);
        const d = await samtools.exec(`view -q 20 ${bam.path} ${locVal}`);
        setData(d);
      }
    })();
  }, [samtools, fileVal, locVal]);

  useEffect(() => {
    const ctx = ref.current.getContext("2d");
    ctx.clearRect(0, 0, width, height);
    if (data) {
      const rows = data?.stdout.split("\n");

      const { start, end } = parseLocString(locVal);
      const bpPerPx = width / (end - start);
      rows.forEach((row, index) => {
        const [name, flagString, refName, startString, mapq, CIGAR] = row.split(
          "\t"
        );
        const start = +startString;
        const flags = +flagString;
        const cigarOps = parseCigar(CIGAR);
        let length = 0;
        for (let i = 0; i < cigarOps.length; i += 2) {
          const len = +cigarOps[i];
          const op = cigarOps[i + 1];
          if (op !== "I" && op !== "S" && op !== "H") {
            length += len;
          }
        }
        if (flags & 2) {
          ctx.fillStyle = "red";
        } else {
          ctx.fillStyle = "blue";
        }
        console.log(start);
        ctx.fillRect(
          +start * bpPerPx,
          index * (featHeight + spacing),
          (start + length) * bpPerPx,
          featHeight
        );
      });
    }
  }, [data, locVal]);
  return (
    <div>
      <p>Enter a BAM file URL:</p>
      <form
        onSubmit={(event) => {
          setFileVal(file);
          setLocVal(loc);
          event.preventDefault();
        }}
      >
        <label htmlFor="url" />
        <input
          id="url"
          type="text"
          value={file}
          onChange={(event) => setFile(event.target.value)}
        />

        <label htmlFor="loc" />
        <input
          id="loc"
          type="text"
          value={loc}
          onChange={(event) => setLoc(event.target.value)}
        />
        <button type="submit">Submit</button>
      </form>
      <canvas ref={ref} width={width} height={height} />
    </div>
  );
}

export default App;
