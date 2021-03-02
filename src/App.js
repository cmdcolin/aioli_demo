import { useCallback, useEffect, useRef, useState } from "react";
import { Aioli } from "@biowasm/aioli";
import GranularRectLayout from "./layout";
import { StringParam, useQueryParams, withDefault } from "use-query-params";

function parseLocString(locString) {
  const [refId, rest] = locString.split(":");
  const [start, end] = rest.split("-");
  return { refId: +refId - 1, start: +start, end: +end };
}

export function parseCigar(cigar) {
  return (cigar || "").split(/([MIDNSHPX=])/);
}

// because we are using use-query-params without a router
export function useForceUpdate() {
  const [, setTick] = useState(0);
  const update = useCallback(() => {
    setTick((tick) => tick + 1);
  }, []);
  return update;
}

const featHeight = 10;
const width = 1800;
const height = 1000;
const initialLoc = "ctgA:1-10000";
const initialFile = "volvox-sorted.bam";
function App() {
  const ref = useRef();
  const [params, setParams] = useQueryParams({
    loc: withDefault(StringParam, initialLoc),
    file: withDefault(StringParam, initialFile),
  });
  const [data, setData] = useState();
  const [samtools, setSamtools] = useState();
  const [bamFile, setBamFile] = useState();
  const [file, setFile] = useState(params.file);
  const [loc, setLoc] = useState(params.loc);
  const forceUpdate = useForceUpdate();

  // this block initializes the samtools tool on Aioli
  useEffect(() => {
    (async () => {
      const tool = new Aioli("samtools/1.10");
      await tool.init();
      setSamtools(tool);
    })();
  }, []);

  // this block "mounts" the BAM file on the Aioli FS
  useEffect(() => {
    (async () => {
      // wait on "samtools" being initialized even though "samtools" is not
      // used here, because it must mount something internally in the
      // Aioli/samtools interface
      if (samtools) {
        const url = new URL(params.file, window.location);
        const bam = await Aioli.mount(`${url}`);
        await Aioli.mount(`${url}.bai`);
        setBamFile(bam);
      }
    })();
  }, [params.file, samtools]);

  // this block performs a query on the bamFile using samtools
  useEffect(() => {
    (async () => {
      if (bamFile && samtools) {
        const d = await samtools.exec(
          `view -q 20 ${bamFile.path} ${params.loc}`
        );
        setData(d);
      }
    })();
  }, [bamFile, params.loc, samtools]);

  // this block draws the rectangles
  useEffect(() => {
    const ctx = ref.current.getContext("2d");
    ctx.clearRect(0, 0, width, height);
    const { start: locStart, end: locEnd } = parseLocString(params.loc);
    const bpPerPx = width / (locEnd - locStart);
    const layout = new GranularRectLayout();
    if (data) {
      const rows = data?.stdout.split("\n");
      console.log({ rows });

      rows.forEach((row, index) => {
        const [, flagString, , startString, , CIGAR] = row.split("\t");
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
        if (flags & 16) {
          ctx.fillStyle = "#f99";
        } else {
          ctx.fillStyle = "#99f";
        }

        const end = start + length;
        const leftPx = (start - locStart) * bpPerPx;
        const width = (end - start) * bpPerPx;
        const y = layout.addRect(index, start, end, featHeight);

        ctx.fillRect(leftPx, y, width, featHeight);
      });
    }
  }, [data, params.loc]);

  return (
    <div>
      <p>Enter a BAM file URL:</p>
      <form
        onSubmit={(event) => {
          setParams({ file, loc });
          setBamFile(undefined);
          forceUpdate();
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
      {!data ? "Loading..." : null}
      <canvas ref={ref} width={width} height={height} />
    </div>
  );
}

export default App;
